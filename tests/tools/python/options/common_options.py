# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import argparse
import logging
import sys
from datetime import datetime

#-------------------------------------------------------------------------------
# Options parsing
#-------------------------------------------------------------------------------

""" Parent parser that contains common options """
def common_options(prog_name = sys.argv[0]):

    parser = argparse.ArgumentParser(add_help=False, prog=prog_name)

    parser.add_argument("--debug", "-d",
                        help=("Enables debug logging"), 
                        action="store_true")
    parser.add_argument("--logfile",
                        help=("Write log in file"), 
                        default = "", type=str)

    return parser


""" Parse common options from commandline """
def parse_common_options(args, prog_name = sys.argv[0]):

    # Set log level
    logLevel = logging.INFO
    if args.debug:
        logLevel = logging.DEBUG

    # Check logfile path
    if args.logfile != "":
        logDirName  = os.path.dirname(args.logfile)

        if logDirName != "" and not os.path.isdir(logDirName):
            os.mkdir(logDirName, 0o755)

        if not os.path.isfile(args.logfile):
            with open(args.logfile, 'x'): pass

        # Write log in logfile
        logging.basicConfig(filename=args.logfile, 
                            format=("%(levelname)s (%(module)s):"
                                    " %(message)s"), level=logLevel)
    else:
        # Write log to terminal
        logging.basicConfig(format=("%(levelname)s (%(module)s):"
                                    " %(message)s"), level=logLevel)

    logging.info("\n--------------------------------------------------------\n"
                 "Running "+ prog_name + " " +\
                 str(datetime.now().strftime("%d/%m/%Y %H:%M"))+\
                 "\n--------------------------------------------------------\n")

    return logging

#-------------------------------------------------------------------------------
