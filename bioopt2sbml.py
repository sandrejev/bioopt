from bioopt_parser import *
import argparse
import libsbml

warnings.simplefilter("ignore")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts model from bioopt format to sbml format')
    parser.add_argument('bioopt', action='store', help='File containing bioopt model')
    parser.add_argument('sbml', action='store', help='Output file prefix')
    parser.add_argument('--compartment', '-c', dest="c_pattern", action='store', help='Regexp to block extract compartment information from metabolite name')
    parser.add_argument('--inf', dest="inf", default=1000, action='store', type=float, help='Infinity value for a new model (default: 1000)')
    parser.add_argument('--level', '-l', dest="level", default=2, action='store', type=int, help='SBML level')
    parser.add_argument('--version', '-v', dest="version", default=4, action='store', type=int, help='SBML version')

    args = parser.parse_args()

    # Read bioopt model
    parser = BiooptParser()
    model = parser.parse_file(args.bioopt)

    sbml = model.sbml(level=args.level, version=args.version, compartment_pattern=args.c_pattern, inf=args.inf)

    libsbml.writeSBMLToFile(sbml, args.sbml)



