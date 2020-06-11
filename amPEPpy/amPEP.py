import sys
import os
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

def score():
    CTD = {'hydrophobicity': {1: ['R', 'K', 'E', 'D', 'Q', 'N'], 2: ['G', 'A', 'S', 'T', 'P', 'H', 'Y'], 3: ['C', 'L', 'V', 'I', 'M', 'F', 'W']},
           'normalized.van.der.waals': {1: ['G', 'A', 'S', 'T', 'P', 'D', 'C'], 2: ['N', 'V', 'E', 'Q', 'I', 'L'], 3: ['M', 'H', 'K', 'F', 'R', 'Y', 'W']},
           'polarity': {1: ['L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'], 2: ['P', 'A', 'T', 'G', 'S'], 3: ['H', 'Q', 'R', 'K', 'N', 'E', 'D']},
           'polarizability': {1: ['G', 'A', 'S', 'D', 'T'], 2: ['C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'], 3: ['K', 'M', 'H', 'F', 'R', 'Y', 'W']},
           'charge': {1: ['K', 'R'], 2: ['A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'], 3: ['D', 'E']},
           'secondary': {1: ['E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'], 2: ['V', 'I', 'Y', 'C', 'W', 'F', 'T'], 3: ['G', 'N', 'P', 'S', 'D']},
           'solvent': {1: ['A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'], 2: ['R', 'K', 'Q', 'E', 'N', 'D'], 3: ['M', 'S', 'P', 'T', 'H', 'Y']}}

if __name__ == "__main__":
    main()
