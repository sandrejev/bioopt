from bioopt_parser import *
from converter import *
import libsbml
import argparse

warnings.simplefilter("ignore")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts model from bioopt format to sbml format')
    parser.add_argument('bioopt', action='store', help='File containing bioopt model')
    parser.add_argument('sbml', action='store', help='Output file prefix')
    parser.add_argument('--compartment', '-c', dest="c_pattern", action='store', help='Regexp to block extract compartment information from metabolite name')
    parser.add_argument('--in-inf', dest="in_inf", default=1000, action='store', type=float, help='Infinity value in parsed model (default: 1000)')
    parser.add_argument('--out-inf', dest="out_inf", default=1000, action='store', type=float, help='Infinity value for a new model (default: 1000)')
    parser.add_argument('--level', '-l', dest="level", default=2, action='store', type=int, help='SBML level')
    parser.add_argument('--version', '-v', dest="version", default=3, action='store', type=int, help='SBML version')
    parser.add_argument('--reaction-id', dest="reaction_id", default="name", action='store', help="Strategy to generate unique reaction id. Specify 'name' to use reaction name as SBML id or 'auto' to use auto-incrementing id R_XXXX. (default: name)")
    parser.add_argument('--metabolite-id', dest="metabolite_id", default="name", action='store', help="Strategy to generate unique metabolite id. Specify 'name' to use metabolite name as SBML id or 'auto' to use auto-incrementing id M_XXXX. (default: name)")
    parser.add_argument('--compartment-id', dest="compartment_id", default="name", action='store', help="Strategy to generate unique compartment id. Specify 'name' to use compartment name as SBML id or 'auto' to use auto-incrementing id C_XXXX. (default: name)")

    args = parser.parse_args()

    # Read bioopt model
    parser = BiooptParser(inf=args.in_inf)
    model = parser.parse_file(args.bioopt)

    converter = Bioopt2SbmlConverter(level=args.level, version=args.version,
          compartment_pattern=args.c_pattern, inf=args.out_inf,
          reaction_id=args.reaction_id, metabolite_id=args.metabolite_id, compartment_id=args.compartment_id)

    sbml = converter.convert(model)
    libsbml.writeSBMLToFile(sbml, args.sbml)

    print "Finished converting {0} into SBML ({1})".format(args.bioopt, args.sbml)



