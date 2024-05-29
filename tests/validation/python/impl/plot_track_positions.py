# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
import plotting

# python includes
import numpy as np
import pandas as pd
import os


""" Read the recorded track positions from files and prepare data frames """
def read_navigation_data(inputdir, read_cuda, logging):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    ray_data_file = ray_data_cuda_file = ""
    helix_data_file = helix_data_cuda_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if read_cuda and filename.find('ray_navigation_track_pos_cuda.') != -1:
            ray_data_cuda_file = inputdir + "/" + filename
        elif filename.find('ray_navigation_track_pos.') != -1:
            ray_data_file = inputdir + "/" + filename
        elif read_cuda and filename.find('helix_navigation_track_pos_cuda.') != -1:
            helix_data_cuda_file = inputdir + "/" + filename
        elif filename.find('helix_navigation_track_pos.') != -1:
            helix_data_file = inputdir + "/" + filename

    # Read navigation data
    def read_data(file):
        if file:
            # Preserve floating point precision
            df = pd.read_csv(file, float_precision='round_trip')
            logging.debug(df)
        else:
            logging.warning("Could not find navigation data file: " + file)
            df = pd.DataFrame({})

        return df

    ray_df = read_data(ray_data_file)
    ray_cuda_df = read_data(ray_data_cuda_file)
    helix_df = read_data(helix_data_file)
    helix_cuda_df = read_data(helix_data_cuda_file)

    return ray_df, ray_cuda_df, helix_df, helix_cuda_df


""" Plot the track positions of two data sources - rz view """
def compare_track_pos_xy(opts, detector, scan_type, plotFactory, out_format,
                         df1, label1, color1, df2, label2, color2):

    n_rays = np.max(df1['track_id']) + 1
    tracks = "rays" if scan_type == "ray" else "helices"

    # Reduce data to the requested z-range (50mm tolerance)
    min_z = opts.z_range[0]
    max_z = opts.z_range[1]
    assert min_z < max_z, "xy plotting range: min z must be smaller that max z"
    pos_range = lambda data: ((data['z'] > min_z) & (data['z'] < max_z))

    first_x, first_y = plotting.filter_data(
                                         data      = df1,
                                         filter    = pos_range,
                                         variables = ['x', 'y'])

    second_x, second_y = plotting.filter_data(
                                         data      = df2,
                                         filter    = pos_range,
                                         variables = ['x', 'y'])

    # Plot the xy coordinates of the filtered track positions
    lgd_ops = plotting.legend_options('upper center', 4, 0.4, 0.005)
    hist_data = plotFactory.scatter(
                            figsize = (8, 8),
                            x       = first_x,
                            y       = first_y,
                            xLabel  = r'$x\,\mathrm{[mm]}$',
                            yLabel  = r'$y\,\mathrm{[mm]}$',
                            label   = label1,
                            color   = color1,
                            alpha   = 1.,
                            showStats = lambda x, y: f"{n_rays} {tracks}",
                            lgd_ops = lgd_ops)

    # Compare agaist second data set
    plotFactory.highlight_region(hist_data, second_x, second_y, color2, label2)

    # Refine legend
    hist_data.lgd.legend_handles[0].set_visible(False)
    for handle in hist_data.lgd.legend_handles[1:]:
        handle.set_sizes([40])

    # For this plot, move the legend ouside
    hist_data.lgd.set_bbox_to_anchor((0.5, 1.11))

    # Adjust spacing in box
    for vpack in hist_data.lgd._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

    detector_name = detector.replace(' ', '_')
    l1 = label1.replace(' ', '_').replace("(", "").replace(")", "")
    l2 = label2.replace(' ', '_').replace("(", "").replace(")", "")

    # Need a very high dpi to reach a good coverage of the individual points
    plotFactory.write_plot(
                        hist_data,
                        f"{detector_name}_{scan_type}_track_pos_{l1}_{l2}_xy",
                        out_format, dpi = 600)


""" Plot the track positions of two data sources - rz view """
def compare_track_pos_rz(opts, detector, scan_type, plotFactory, out_format,
                         df1, label1, color1, df2, label2, color2):

    n_rays = np.max(df1['track_id']) + 1
    tracks =  "rays" if scan_type == "ray" else "helices"

    first_x, first_y, first_z = plotting.filter_data(
                                                data      = df1,
                                                variables = ['x', 'y', 'z'])

    second_x, second_y, second_z = plotting.filter_data(
                                                data      = df2,
                                                variables = ['x', 'y', 'z'])

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = plotting.legend_options('upper center', 4, 0.8, 0.005)
    hist_data = plotFactory.scatter(
                            figsize = (12, 6),
                            x      = first_z,
                            y      = np.hypot(first_x, first_y),
                            xLabel = r'$z\,\mathrm{[mm]}$',
                            yLabel = r'$r\,\mathrm{[mm]}$',
                            label  = label1,
                            color  = color1,
                            alpha   = 1.,
                            showStats = lambda x, y: f"{n_rays} {tracks}",
                            lgd_ops = lgd_ops)

    # Compare agaist second data set
    plotFactory.highlight_region(hist_data, second_z, np.hypot(second_x, second_y), color2, label2)

    # Refine legend
    hist_data.lgd.legend_handles[0].set_visible(False)
    for handle in hist_data.lgd.legend_handles[1:]:
        handle.set_sizes([40])

    # For this plot, move the legend ouside
    hist_data.lgd.set_bbox_to_anchor((0.5, 1.15))

    # Adjust spacing in box
    for vpack in hist_data.lgd._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

    detector_name = detector.replace(' ', '_')
    l1 = label1.replace(' ', '_').replace("(", "").replace(")", "")
    l2 = label2.replace(' ', '_').replace("(", "").replace(")", "")

    # Need a very high dpi to reach a good coverage of the individual points
    plotFactory.write_plot(
                        hist_data,
                        f"{detector_name}_{scan_type}_track_pos_{l1}_{l2}_rz",
                        out_format, dpi = 600)


""" Plot the absolute track positions distance """
def plot_track_pos_dist(opts, detector, scan_type, plotFactory, out_format,
                        df1, label1, df2, label2):

    n_rays = np.max(df1['track_id']) + 1
    tracks =  "rays" if scan_type == "ray" else "helices"

    dist = np.sqrt(np.square(df1['x'] - df2['x']) +
                   np.square(df1['y'] - df2['y']) +
                   np.square(df1['z'] - df2['z']))

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    hist_data = plotFactory.hist1D(x       = dist,
                                   bins    = 100,
                                   xLabel  = r'$d\,\mathrm{[mm]}$',
                                   setLog  = True,
                                   lgd_ops = lgd_ops)

    # Adjust spacing in box
    hist_data.lgd.legend_handles[0].set_visible(False)
    for vpack in hist_data.lgd._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

    detector_name = detector.replace(' ', '_')
    l1 = label1.replace(' ', '_').replace("(", "").replace(")", "")
    l2 = label2.replace(' ', '_').replace("(", "").replace(")", "")
    plotFactory.write_plot(hist_data,
                           f"{detector_name}_{scan_type}_dist_{l1}_{l2}",
                           out_format)


""" Plot the track position residual for the given variable """
def plot_track_pos_res(opts, detector, scan_type, plotFactory, out_format,
                       df1, label1, df2, label2, var):

    n_rays = np.max(df1['track_id']) + 1
    tracks =  "rays" if scan_type == "ray" else "helices"

    res = df1[var] - df2[var]

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    hist_data = plotFactory.hist1D(x       = res,
                                   bins    = 100,
                                   xLabel  = r'$\mathrm{res}' + rf'\,{var}' + r'\,\mathrm{[mm]}$',
                                   setLog  = True,
                                   lgd_ops = lgd_ops)

    # Adjust spacing in box
    hist_data.lgd.legend_handles[0].set_visible(False)
    for vpack in hist_data.lgd._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

    detector_name = detector.replace(' ', '_')
    l1 = label1.replace(' ', '_').replace("(", "").replace(")", "")
    l2 = label2.replace(' ', '_').replace("(", "").replace(")", "")
    plotFactory.write_plot(hist_data,
                           f"{detector_name}_{scan_type}_res_{var}_{l1}_{l2}",
                           out_format)
