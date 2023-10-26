# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import plot_ray_scan
from pyplot_factory import pyplot_factory

# python includes
import argparse
import logging
import numpy as np
import pandas as pd
import os
import sys
from datetime import datetime
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def __main__():

#----------------------------------------------------------------arg parsing

    parser = argparse.ArgumentParser(description = "Detray Ray Scan")
    parser.add_argument("--debug", "-d",
                        help=("Enables debug logging"), 
                        action="store_true")
    parser.add_argument("--logfile",
                        help=("Write log in file"), 
                        default = "", type=str)
    parser.add_argument("--input", "-i",
                        help=("Input ray scan data file."),
                        default = "", type=str, required=True)
    parser.add_argument("--outdir", "-o",
                        help=("Output directory for plots."),
                        default = "./geometry_plots/", type=str)
    parser.add_argument("--output_format", "-of",
                        help=("Format of the plot files (svg|png|pdf)."),
                        default = "png", type=str)
    parser.add_argument("--z_range", "-zrng", nargs=2,
                        help=("z range for the xy-view."),
                        default = [-50, 50], type=float)
    parser.add_argument("--hide_portals",
                        action="store_true", default=False)
    parser.add_argument("--hide_passives",
                        action="store_true", default=False)

    args = parser.parse_args()

#---------------------------------------------------------------------config

    # Check output path
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir, 0o755)
    outdir = args.outdir

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
                 "Running ray scan validation "+\
                 str(datetime.now().strftime("%d/%m/%Y %H:%M"))+\
                 "\n--------------------------------------------------------\n")

    # Check input data files from material scan
    if args.input == "":
        logging.error(f"Please specify an input data file!")
        sys.exit(1)

    if not os.path.isfile(args.input):
        logging.error(f"Data file does not exist! ({args.input})")
        sys.exit(1)

    if not args.output_format in ["svg", "png", "pdf"]:
        logging.error(f"Unknown output file format: {out_format}")
        sys.exit(1)

    ray_scan_file = args.input
    out_format = args.output_format

#----------------------------------------------------------------prepare data

    # Get detector name
    if ray_scan_file.find('ray_scan_') != -1:
        detector_name = ray_scan_file.removeprefix('ray_scan_')
        scan_type = "ray"
    elif ray_scan_file.find('helix_scan_') != -1:
        detector_name = ray_scan_file.removeprefix('helix_scan_')
        scan_type = "helix"
    else:
        logging.error('Input filename needs to contain \'ray_scan\' prefix')
    detector_name = detector_name.removesuffix('.csv')
    detector_name = detector_name.replace('_', ' ')

    df = pd.read_csv(ray_scan_file)

    plot_factory = pyplot_factory(outdir + "geometry_", logging)

#------------------------------------------------------------------------run

    plot_ray_scan.intersection_points_xy(args, df, detector_name,
                                         scan_type, plot_factory, out_format)
    plot_ray_scan.intersection_points_rz(args, df, detector_name, scan_type,
                                         plot_factory, out_format)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

#------------------------------------------------------------------------------- 
