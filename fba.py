import cplex_utils
import gtrass
from bioopt_parser import BiooptParser
from cplex import infinity as cplex_inf
import itertools
import argparse
from operator import attrgetter as _a
from operator import itemgetter as _i
import re
import sys


def chunkify(lst, n, jumpy=False):
    if jumpy:
        return [lst[i::n] for i in xrange(n)]

    avg = len(lst) / float(n)
    out = []
    last = 0.0

    while last < len(lst):
        out.append(lst[int(last):int(last + avg)])
        last += avg

def geneko2reactionko(knockouts_g, associations, genes):
        genes_values = {g: 1 for g in genes}

        ga_cache = {rxn: {} for rxn in associations}
        for ko_i, ko_genes in enumerate(knockouts_g):
            for g in ko_genes: genes_values[g] = 0
            ko_reactions = set()
            for rxn, ga in associations.iteritems():
                ga_cache_rxn = sum((g in ko_genes)^i for i, g in enumerate(ga.genes))
                if ga_cache_rxn in ga_cache[rxn]:
                    r_active = ga_cache[rxn][ga_cache_rxn]
                else:
                    r_active = ga.eval(genes_values)
                    ga_cache[rxn][ga_cache_rxn] = r_active

                if not r_active:
                    ko_reactions.add(cplex.rxn2i[rxn])

            yield ko_reactions
            for g in ko_genes: genes_values[g] = 1

def main_knockouts(cplex, args):
    # Find reactions candidates for knockouts
    excluded_ko_reactions = set(r.name for r in bioopt.find_reactions("(xtI|xtO|xtX)", True))
    reactions = sorted([r for r in bioopt.reactions if r.name not in excluded_ko_reactions], key=_a('name'))
    reactions_i = [cplex.rxn2i[r.name] for r in reactions]

    if args.gtrass:
        reactions, genes, associations = gtrass.parse_file(args.gtrass)
        if re.match("^\d+$", args.knockouts):
            knockouts_g = itertools.combinations(genes, int(args.knockouts))
        else:
            knockouts_g = [set(l.strip().split("\t")) for l in open(args.knockouts)]

        knockouts_i = geneko2reactionko(knockouts_g, associations, genes)
    else:
        if re.match("^\d+$", args.knockouts):
            knockouts_i = itertools.combinations(reactions_i, int(args.knockouts))
        else:
            knockouts_i = list()
            for l in open(args.knockouts):
                knockouts_i.append(set(l.strip().split("\t")))

        knockouts_g = [""]*len(knockouts_i)

    # Run knockouts
    print "Genes\tReactions\tStatus\tObjective"
    for ko_genes, ko_res in itertools.izip(knockouts_g, loop_knockouts(cplex, knockouts_i)):
        ko_i, ko_status, ko_obj = ko_res
        ko_rxn = ",".join([cplex.i2rxn[r_i] for r_i in ko_i])
        ko_genes = ",".join(ko_genes)
        print "{}\t{}\t{}\t{}".format(ko_genes, ko_rxn, ko_status, ko_obj)


def loop_knockouts(cplex, knockouts_i):
    cm = cplex.model
    cm.parameters.lpmethod.set(cm.parameters.lpmethod.values.primal)
    cm.parameters.preprocessing.reduce.set(3)

    lb_default = {r_i: -cplex_inf if b.lb <= -cplex_inf else b.lb for r_i, b in cplex.i2bounds.iteritems()}
    ub_default = {r_i: cplex_inf if b.ub >= cplex_inf else b.ub for r_i, b in cplex.i2bounds.iteritems()}

    for i, rxns in enumerate(knockouts_i):
        new_bounds, bck_ub, bck_lb = [], [], []
        for r_i in rxns:
            cplex_bounds = (r_i, 0.0)
            new_bounds.append(cplex_bounds)
            bck_lb.append((r_i, lb_default[r_i]))
            bck_ub.append((r_i, ub_default[r_i]))

        if len(rxns):
            cm.variables.set_lower_bounds([(r_i, 0.0) for r_i in rxns])
            cm.variables.set_upper_bounds([(r_i, 0.0) for r_i in rxns])

        cm.solve()
        status = cm.solution.get_status_string()
        value = cm.solution.get_objective_value() if cplex_utils.is_optimal(cm) else 0.0

        yield rxns, status, value

        if len(rxns):
            cm.variables.set_lower_bounds(bck_lb)
            cm.variables.set_upper_bounds(bck_ub)

def main_simple(cplex, args):
    cm = cplex.model
    cm.set_log_stream(sys.stdout)
    cm.set_error_stream(sys.stdout)
    cm.set_warning_stream(sys.stdout)
    cm.set_results_stream(sys.stdout)
    cm.solve()

    status = cm.solution.get_status_string()
    value = cm.solution.get_objective_value() if cplex_utils.is_optimal(cm) else 0.0

    print "Solution: {} ({})".format(value, status)

    if cplex_utils.is_optimal(cm):
        if args.show_fluxes:
            show_all_fluxes = args.show_fluxes == "*"
            if not show_all_fluxes:
                re_show_fluxes = re.compile(args.show_fluxes, re.M | re.I)

            primal = sorted([(cplex.i2rxn[r_i], val) for r_i, val in enumerate(cm.solution.get_values())], key=_i(0))

            print "Fluxes\n==============================="
            for rxn, val in primal:
                if val != 0 and (show_all_fluxes or re_show_fluxes.search(rxn)):
                    print "{:<20}\t{}".format(rxn, val)
            print "==============================="

        if args.show_dual:
            print "Dual values\n==============================="
            show_all_dual = args.show_dual == "*"
            if not show_all_dual:
                re_show_dual = re.compile(args.show_dual, re.M | re.I)

            dual = sorted([(cplex.i2cpd[c_i], val) for c_i, val in enumerate(cm.solution.get_dual_values())], key=_i(0))
            print "Dual\n==============================="
            for cpd, val in dual:
                if val != 0 and (show_all_dual or re_show_dual.search(cpd)):
                    print "{:<20}\t{}".format(cpd, val)
            print "==============================="


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('model', action='store', help='')
    parser.add_argument('objective', action='store', help='')
    parser.add_argument('--gtrass', dest='gtrass', required=False, action='store', help="Gene<->Reaction association file", default=None)
    parser.add_argument('--knockouts', dest='knockouts', required=False, action='store', help="Either file or number of knockouts", default=0)
    parser.add_argument('--show-fluxes', dest='show_fluxes', required=False, action='store')
    parser.add_argument('--show-dual', dest='show_dual', required=False, action='store')
    args = parser.parse_args()

    # Read model
    parser = BiooptParser()
    bioopt = parser.parse_file(args.model)
    if not bioopt.find_reaction(args.objective):
        print "Objective {} not found".format(args.objective)
        exit()

    cplex = cplex_utils.bioopt2cplex(bioopt, objective=args.objective)
    cplex.model.objective.set_linear(cplex.rxn2i[args.objective], 1.0)

    if args.knockouts:
        main_knockouts(cplex, args)
    else:
        main_simple(cplex, args)

