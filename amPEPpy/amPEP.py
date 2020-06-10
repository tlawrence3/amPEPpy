import sys
import argparse
import numpy
from Bio import SeqIO
from amPEPpy._version import __version__

def main():
    parser = argparse.ArgumentParser(prog="ampep", add_help=False)
    parser.add_argument('-v', '--version', action='version', version="v{}".format(__version__))
    parser.add_argument("--seed", type=int, default=os.urandom(64), dest="seed",
                        help="Seed number as integer. This allows reproducibility of output community. Default to random seed number.")
    subparsers = parser.add_subparsers(help="sub-command help")

    parser_train = subparsers.add_parser("train", parents=[parser], help='')
    parser_train.set_defaults(func=train)

    parser_classify = subparsers.add_parser("classify", parents=[parser], help="")
    parser_classify.set_defaults(func=classify)

    if len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        args = parser.parse_args()
        args.func(args)

def train(args):
    pass

def classify(args):
    pass

if __name__ == "__main__":
    main()
