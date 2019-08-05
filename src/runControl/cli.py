# -*- coding: utf-8 -*-
"""The command line file

"""

import argparse
import sys
from pkg_resources import (get_distribution, DistributionNotFound)

try:
    __version__ = get_distribution('runControl').version
except DistributionNotFound:
    __version__ = 'unknown'

# Parse command line
parser = argparse.ArgumentParser()
# First define all option groups
group1 = parser.add_argument_group('Input file', 'Required input')
group1.add_argument("-f", "--csv", default="", type=str, dest="f",
                    help="input reads in csv format")
group1.add_argument('-v', '--version', action='version',
                    version=__version__)

# Exit if no input file is provided
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()


def main(args=None):

    import logging
    import logging.handlers

    args = parser.parse_args(args=args)

    log_format = '%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s'
    logging.basicConfig(filename='runko.log', level=logging.INFO,
                        format=log_format, datefmt='%Y/%m/%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    from runControl import run
    run.main(args.f)
