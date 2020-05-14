#! /usr/bin/env python

from .parser import create_parser
from .run import ebfilter_main

def main():

    parser = create_parser()
    args = parser.parse_args()
    ebfilter_main(args)

