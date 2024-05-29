# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
from impl import plot_material_scan as mat_plotter
from plotting import pyplot_factory as plt_factory
from options import common_options, plotting_options
from options import parse_common_options, parse_plotting_options

# python includes
import argparse
import pandas as pd


def __main__():

#----------------------------------------------------------------arg parsing

    descr = "Detray Material Validation"

    common_parser = common_options(descr)
    plotting_parser = plotting_options()

    parser = argparse.ArgumentParser(description = descr, parents=[common_parser, plotting_parser])

    args = parser.parse_args()

    logging = parse_common_options(args, descr)
    input_dir, out_dir, out_format = parse_plotting_options(args, logging)

#----------------------------------------------------------------prepare data

    # Get detector name
    detector_name = mat_scan_file.removeprefix('material_scan_')
    detector_name = detector_name.removesuffix('.csv')
    detector_name = detector_name.replace('_', ' ')

    df = pd.read_csv(mat_scan_file)

    plot_factory = plt_factory(outdir + "material_", logging)

#------------------------------------------------------------------------run

    # The histograms are not re-weighted (if the rays are not evenly distributed
    # the material in some bins might be artificially high)!
    mat_plotter.X0_vs_eta_phi(df, detector_name, plot_factory, out_format)
    mat_plotter.L0_vs_eta_phi(df, detector_name, plot_factory, out_format)
    mat_plotter.X0_vs_eta(df, detector_name, plot_factory, out_format)
    mat_plotter.L0_vs_eta(df, detector_name, plot_factory, out_format)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

#------------------------------------------------------------------------------- 
