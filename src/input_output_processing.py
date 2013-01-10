#!/usr/bin/env python
import fileinput
import sys
from argparse import ArgumentParser

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
        help='specify input files, default is stdin')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-o', '--output', 
        help='specify the output file.  The default is stdout')
    args = parser.parse_args()

    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')

    for line in fileinput.input(args.files, inplace=args.inplace):
        print line,
