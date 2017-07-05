from bioopt_parser import *
import argparse
import numpy as np


def bioopt2optflux(model):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts moel from bioopt format to optflux format')
    parser.add_argument('bioopt', action='store', help='File containing bioopt model')
    parser.add_argument('optflux', action='store', help='Output file prefix')

    args = parser.parse_args()

    # Read bioopt model
    parser = BiooptParser()
    model = parser.parse_file(args.bioopt)

    # Write list of metabolites
    metabolites = model.find_metabolites(None)
    f_mets = open(args.optflux + '.mets', 'w')
    for m in metabolites:
        f_mets.write("{0}\n".format(m.name))
    f_mets.close()

    # Write list of constraints
    f_constr = open(args.optflux + '.constr', 'w')
    for r in model.reactions:
        b = r.find_effective_bounds()
        f_constr.write("{0}\t{1}\t{2}\n".format(r.name, b.lb, b.ub))
    f_constr.close()

    # Write stoichiometry matrix
    # TODO:

    stoich = np.zeros((len(metabolites), len(model.reactions)))
    for r_i, r in enumerate(model.reactions):
        r_metabolites = [mb.metabolite for mb in r.reactants]
        p_metabolites = [mb.metabolite for mb in r.products]

        for m_i, m in enumerate(metabolites):
            #if m_i > 0:
            #    f_matrix.write("\t")

            coef = 0
            if m in r_metabolites:
                coef = -r.reactants.find_member(m.name).coefficient
            elif m in p_metabolites:
                coef = r.products.find_member(m.name).coefficient

            stoich[m_i, r_i] = coef


    with open(args.optflux + '.stoich', 'w') as f_matrix:
        f_matrix.writelines('\t'.join([str(j) for j in i]) + '\n' for i in stoich)
    f_matrix.close()