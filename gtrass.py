from __future__ import with_statement
import ast
import re
from sympy import SOPform, POSform, And, Or
import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def _time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException, "Timed out!"
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

def _eval(association_tree, variables):
    if isinstance(association_tree, ast.Num): # <number>
        return association_tree.n

    if isinstance(association_tree, ast.Name):
        if association_tree.id == "None":
            return None
        else:
            return variables[association_tree.id]

    elif isinstance(association_tree, ast.BoolOp):
        if isinstance(association_tree.op, ast.And):
            return all(_eval(v, variables) for v in association_tree.values)
        if isinstance(association_tree.op, ast.Or):
            return any(_eval(v, variables) for v in association_tree.values)
    else:
        raise TypeError(association_tree)


def _sympy_str(exp, first=True):
    brackets = "{}" if first else "({})"
    if isinstance(exp, Or) or isinstance(exp, And):
        if isinstance(exp, Or):
            op = " or "
        if isinstance(exp, And):
            op = " and "
        return brackets.format(op.join(_sympy_str(arg, first=False) for arg in exp.args))
    else:
        return exp.name

re_gene_sep = re.compile(r"(?:\(|\)| and | or |;| )")
re_or = re.compile(" or ")
re_and = re.compile(" and ")

class GeneAssociation():
    def __init__(self, association_str, association_tree):
        self.genes = set(g for g in re_gene_sep.split(association_str) if g)
        self.association_str = association_str
        self.association_tree = association_tree

    def eval(self, variables):
        return _eval(self.association_tree, variables)

    def simplify(self, name=""):
        match_or = re_or.search(self.association_str)
        match_and = re_and.search(self.association_str)
        genes = list(self.genes)

        if len(genes) == 1:
            return parse(genes[0])
        elif match_or and not match_and:
            return parse(" or ".join(genes))
        elif match_and and not match_or:
            return parse(" and ".join(genes))
        else:
            genes_mat_size = len(genes)

            pos_timeout, sop_timeout = True, True
            if len(genes) < 10:
                genes_mat = [[int(x) for x in '{:0b}'.format(x).zfill(genes_mat_size)] for x in range(2**genes_mat_size)]
                genes_mat = [values for values in genes_mat if self.eval(dict(zip(genes, values)))]

                try:
                    with _time_limit(10):
                        association_pos = _sympy_str(POSform(genes, genes_mat))
                        pos_timeout = False
                except TimeoutException:
                    association_pos = self.association_str

                try:
                    with _time_limit(10):
                        association_sop = _sympy_str(SOPform(genes, genes_mat))
                        sop_timeout = False
                except TimeoutException:
                    association_sop = self.association_str
            else:
                association_pos = self.association_str
                association_sop = self.association_str

            ga_str = [self.association_str, association_pos, association_sop]
            ga_tree = parse(min(ga_str, key=len))

            if sop_timeout and pos_timeout:
                print "{}\t{}".format(name, self.association_str)

            return ga_tree



def parse(association_str):
    association_str = re.sub(" *; *", " or ", association_str)
    tree = ast.parse(association_str, "", "eval").body
    return GeneAssociation(association_str, tree)

def _map(tree, gene_map={}):
    if isinstance(tree, ast.Name):
        return gene_map.get(tree.id, tree.id)
    elif isinstance(tree, ast.BoolOp):
        op = " or " if isinstance(tree.op, ast.Or) else " and "
        return op.join([_map(v) for v in tree.values])

def map(association_str, gene_map={}):
    tree = ast.parse(association_str, "", "eval").body
    return _map(tree, gene_map)


def parse_file(path):
    genes, reactions, associations = set(), set(), {}
    for l in open(path):
        rxn, ga = l.strip().split("\t")
        associations[rxn] = parse(ga)
        reactions.add(rxn)
        for g in re_gene_sep.split(ga):
            genes.add(g)

    return sorted(reactions), sorted(genes), associations
