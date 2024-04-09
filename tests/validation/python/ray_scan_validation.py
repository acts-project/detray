# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
from validation import plot_ray_scan as scan_plotter
from validation import plt_factory

# python includes
import argparse
import logging
import numpy as np
import pandas as pd
import os
import sys
from datetime import datetime


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
                        help=("Hide intersections with portal surfaces."),
                        action="store_true", default=False)
    parser.add_argument("--hide_passives",
                        help=("Hide intersections with passive surfaces."),
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
    if ray_scan_file.find('_ray_scan') != -1:
        detector_name = ray_scan_file.removesuffix('_ray_scan.csv')
        scan_type = "ray"
    elif ray_scan_file.find('_helix_scan') != -1:
        detector_name = ray_scan_file.removesuffix('_helix_scan.csv')
        scan_type = "helix"
    else:
        logging.error('Input filename needs to contain \'ray_scan\' suffix')
    detector_name = detector_name.replace('_', ' ')

    df = pd.read_csv(ray_scan_file)

    plot_factory = plt_factory(outdir + "geometry_", logging)

#------------------------------------------------------------------------run

    scan_plotter.intersection_points_xy(args, df, detector_name,
                                        scan_type, plot_factory, out_format)
    scan_plotter.intersection_points_rz(args, df, detector_name, scan_type,
                                        plot_factory, out_format)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

#------------------------------------------------------------------------------- 
