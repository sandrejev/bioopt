import ast
import re

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


re_gene_sep = re.compile(r"(?:\(|\)| and | or |;| )")

class GeneAssociation():
    def __init__(self, association_str, association_tree):
        self.genes = set(re_gene_sep.split(association_str))
        self.association_str = association_str
        self.association_tree = association_tree

    def eval(self, variables):
        return _eval(self.association_tree, variables)


def parse(association_str):
    association_str = re.sub(" *; *", " or ", association_str)
    tree = ast.parse(association_str, "", "eval").body
    return GeneAssociation(association_str, tree)

def parse_file(path):
    genes, reactions, associations = set(), set(), {}
    for l in open(path):
        rxn, ga = l.strip().split("\t")
        associations[rxn] = parse(ga)
        reactions.add(rxn)
        for g in re_gene_sep.split(ga):
            genes.add(g)

    return sorted(reactions), sorted(genes), associations
