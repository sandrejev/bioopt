from bioopt_parser import *
import argparse
import os
import numpy as np


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts moel from bioopt format to optflux format')
    parser.add_argument('bioopt', action='store', nargs='+', help='Files containing bioopt models')
    parser.add_argument('output', action='store', help='Output file')

    args = parser.parse_args()

    models = []
    for bioopt_path in args.bioopt:
        print "Parsing {0}: ".format(bioopt_path),
        if not os.path.exists(bioopt_path):
            print "File not found!"
            break

        parser = BiooptParser()
        model = parser.parse_file(bioopt_path)
        models.append(model)
        print "Done"

    com_model = Model.commune(models)
    print "Community model created"

    com_model.save(args.output)
    print "File {0} written!".format(args.output)