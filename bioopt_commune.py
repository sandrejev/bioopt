from bioopt_parser import *
import argparse
import os
import numpy as np


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts moel from bioopt format to optflux format')
    parser.add_argument('bioopt', action='store', nargs='+', help='Files containing bioopt models')
    parser.add_argument('output', action='store', help='Output file')
    parser.add_argument('--block', '-b', dest="block", action='store', help='Regexp to block uptake of reactions')
    parser.add_argument('--inf', dest="inf", default=None, action='store', help='Infinity value for a new model (default: 1000)')

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


    if args.inf is None:
        args.inf = 1000

    com_model = Model.commune(models)
    print "Community model created"

    if args.block:
        print "Blocking reactions: ",
        block_re = re.compile(args.block)
        blocked_rcount = 0
        for r in com_model.reactions:
            if block_re.match(r.name):
                if blocked_rcount > 0:
                    print ", ",
                print r.name,
                r.bounds = Bounds(0, 0)
                blocked_rcount += 1
        print ""

    for r in model.reactions:
        if math.fabs(r.bounds.lb) == Bounds.inf():
            r.bounds.lb = args.inf if r.bounds.lb >= 0 else -args.inf

        if math.fabs(r.bounds.ub) >= Bounds.inf():
            r.bounds.ub = args.inf if r.bounds.ub >= 0 else -args.inf

    com_model.save(args.output)
    print "File {0} written!".format(args.output)