# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from .plot_ray_scan import intersection_points_xy, intersection_points_rz

# python includes
import pandas as pd
import os


""" Read the navigation data from files and prepare data frames """
def read_navigation_data(inputdir, logging):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    detector_name = "default_detector"
    ray_scan_intersections_file = ray_scan_track_param_file = ""
    helix_scan_intersections_file = helix_scan_track_param_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if filename.find('_ray_scan_intersections') != -1:
            detector_name = ray_scan_intersections_file.removesuffix('_ray_scan_intersections.csv')
            ray_scan_intersections_file = inputdir + "/" + filename
        elif filename.find('_ray_scan_track_parameters') != -1:
            ray_scan_track_param_file = inputdir + "/" + filename
        elif filename.find('_helix_scan_intersections') != -1:
            detector_name = helix_scan_intersections_file.removesuffix('_helix_scan_intersections.csv')
            helix_scan_intersections_file = inputdir + "/" + filename
        elif filename.find('_helix_scan_track_parameters') != -1:
            helix_scan_track_param_file = inputdir + "/" + filename

    detector_name = detector_name.replace('_', ' ')

    # Read ray scan data
    if ray_scan_intersections_file:
        ray_inters_df = pd.read_csv(ray_scan_intersections_file)
        ray_trk_param_df = pd.read_csv(ray_scan_track_param_file)
        ray_scan_df = pd.concat([ray_inters_df, ray_trk_param_df], axis=1)

        logging.debug(ray_scan_df)
    else:
        logging.warning("Could not find ray scan data: " + ray_scan_intersections_file)
        ray_scan_df = pd.DataFrame({})

    # Read helix scan data
    if helix_scan_intersections_file:
        helix_inters_df = pd.read_csv(helix_scan_intersections_file)
        helix_trk_param_df = pd.read_csv(helix_scan_track_param_file)
        helix_scan_df = pd.concat([helix_inters_df, helix_trk_param_df], axis=1)

        logging.debug(helix_scan_df)
    else:
        logging.warning("Could not find helix scan data" + helix_scan_intersections_file)
        helix_scan_df = pd.DataFrame({})

    return detector_name, ray_scan_df, helix_scan_df



""" Plot the data gathered during the navigaiton validation """
def plot_navigation_data(args, detector_name, plot_factory, df_truth, data_name, out_format = "png"):

    # Plot truth scan
    intersection_points_xy(args, df_truth, detector_name,
                                        data_name, plot_factory, out_format)
    intersection_points_rz(args, df_truth, detector_name,
                                        data_name, plot_factory, out_format)
