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
import math

""" Read track position data """
def read_track_data(file, logging):
    if file:
        # Preserve floating point precision
        df = pd.read_csv(file, float_precision='round_trip')
        logging.debug(df)
    else:
        logging.warning("Could not find navigation data file: " + file)
        df = pd.DataFrame({})

    return df


""" Plot the distributions of track parameter data """
def plot_track_params(opts, detector, track_type, plotFactory, out_format,
                      df):

    from matplotlib.ticker import ScalarFormatter

    detector_name = detector.replace(' ', '_')
    n_tracks = len(df['track_id'])
    tracks =  "rays" if track_type == "ray" else "helices"
    # Where to place the legend box
    box_anchor_x = 1.02
    box_anchor_y = 1.245

    # Plot the charge
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    q_hist_data = plotFactory.hist1D(x       = df['q'],
                                     bins    = 4,
                                     xLabel  = r'$q\,\mathrm{[e]}$',
                                     lgd_ops = lgd_ops,
                                     figsize  = (8.5, 8.5),
                                     ax_formatter = ScalarFormatter())

    # For this plot, move the legend ouside
    q_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(q_hist_data,
                           f"{detector_name}_{track_type}_charge_dist",
                           out_format)

    # Plot the total momentum
    p = np.sqrt(np.square(df['px']) +  np.square(df['py']) + 
                np.square(df['pz']))

    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    p_hist_data = plotFactory.hist1D(x       = p,
                                     bins    = 100,
                                     xLabel  = r'$p_{tot}\,\mathrm{[GeV]}$',
                                     lgd_ops = lgd_ops,
                                     setLog  = True,
                                     figsize  = (8.5, 8.5),
                                     ax_formatter = ScalarFormatter())

    p_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(p_hist_data,
                           f"{detector_name}_{track_type}_p_dist",
                           out_format)

    # Plot the transverse momentum
    pT = np.sqrt(np.square(df['px']) + np.square(df['py']))

    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    pT_hist_data = plotFactory.hist1D(x      = pT,
                                      bins    = 100,
                                      xLabel  = r'$p_{T}\,\mathrm{[GeV]}$',
                                      lgd_ops = lgd_ops,
                                      figsize  = (8.5, 8.5),
                                      ax_formatter = ScalarFormatter())

    pT_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(pT_hist_data,
                           f"{detector_name}_{track_type}_pT_dist",
                           out_format)

    # Plot the x-origin
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    x_hist_data = plotFactory.hist1D(x      = df['x'],
                                     bins    = 100,
                                     xLabel  = r'$x\,\mathrm{[mm]}$',
                                     lgd_ops = lgd_ops,
                                     figsize  = (8.5, 8.5),
                                     ax_formatter = ScalarFormatter())

    x_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(x_hist_data,
                           f"{detector_name}_{track_type}_x_origin",
                           out_format)

    # Plot the y-origin
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    y_hist_data = plotFactory.hist1D(x      = df['y'],
                                     bins    = 100,
                                     xLabel  = r'$y\,\mathrm{[mm]}$',
                                     lgd_ops = lgd_ops,
                                     figsize  = (8.5, 8.5),
                                     ax_formatter = ScalarFormatter())

    y_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(y_hist_data,
                           f"{detector_name}_{track_type}_y_origin",
                           out_format)
    # Plot the z-origin
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    z_hist_data = plotFactory.hist1D(x      = df['z'],
                                     bins    = 100,
                                     xLabel  = r'$z\,\mathrm{[mm]}$',
                                     lgd_ops = lgd_ops,
                                     figsize  = (8.5, 8.5),
                                     ax_formatter = ScalarFormatter())

    z_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(z_hist_data,
                           f"{detector_name}_{track_type}_z_origin",
                           out_format)

    # Plot the phi angle of the track direction
    phi = np.arctan2(df['py'], df['px'])
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    dir_phi_hist_data = plotFactory.hist1D(x      = phi,
                                        bins    = 100,
                                        xLabel  = r'$\varphi\,\mathrm{[rad]}$',
                                        lgd_ops = lgd_ops,
                                        figsize  = (8.5, 8.5),
                                        ax_formatter = ScalarFormatter())

    dir_phi_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(dir_phi_hist_data,
                           f"{detector_name}_{track_type}_dir_phi",
                           out_format)

    # Plot the theta value of the track direction
    theta = np.arctan2(pT, df['pz'])
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    dir_theta_hist_data = plotFactory.hist1D(x      = theta,
                                             bins    = 100,
                                             xLabel  = r'$\theta\,\mathrm{[rad]}$',
                                             lgd_ops = lgd_ops,
                                             figsize  = (8.5, 8.5),
                                             ax_formatter = ScalarFormatter())

    dir_theta_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(dir_theta_hist_data,
                           f"{detector_name}_{track_type}_dir_theta",
                           out_format)

    # Plot the eta value of the track direction
    eta = np.arctanh(df['pz'] / p)
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    dir_eta_hist_data = plotFactory.hist1D(x       = eta,
                                           bins    = 100,
                                           xLabel  = r'$\eta$',
                                           lgd_ops = lgd_ops,
                                           figsize  = (8.5, 8.5),
                                           ax_formatter = ScalarFormatter())

    dir_eta_hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(dir_eta_hist_data,
                           f"{detector_name}_{track_type}_dir_eta",
                           out_format)


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
                            figsize = (10, 10),
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

    # For this plot, move the legend ouside
    hist_data.lgd.set_bbox_to_anchor((0.5, 1.095))

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

    # For this plot, move the legend ouside
    hist_data.lgd.set_bbox_to_anchor((0.5, 1.168))

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

    n_tracks = np.max(df1['track_id']) + 1
    tracks =  "rays" if scan_type == "ray" else "helices"
    # Where to place the legend box
    box_anchor_x = 1.02
    box_anchor_y = 1.24

    dist = np.sqrt(np.square(df1['x'] - df2['x']) +
                   np.square(df1['y'] - df2['y']) +
                   np.square(df1['z'] - df2['z']))

    dist_outlier = math.sqrt(3 * opts.outlier**2)

    # Remove outliers
    filter_dist = np.absolute(dist) < dist_outlier 
    filtered_dist = dist[filter_dist]

    if not np.all(filter_dist == True):
        print(f"\nRemoved outliers (dist):")
        for i, d in enumerate(dist):
            if math.fabs(d) > dist_outlier:
                track_id = (df1['track_id'].to_numpy())[i]
                print(f"track {track_id}: {d}")

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = plotting.legend_options('upper right', 4, 0.8, 0.005)
    hist_data = plotFactory.hist1D(x       = filtered_dist,
                                   bins    = 100,
                                   xLabel  = r'$d\,\mathrm{[mm]}$',
                                   setLog  = True,
                                   figsize  = (8.5, 8.5),
                                   lgd_ops = lgd_ops)

    hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    detector_name = detector.replace(' ', '_')
    l1 = label1.replace(' ', '_').replace("(", "").replace(")", "")
    l2 = label2.replace(' ', '_').replace("(", "").replace(")", "")
    plotFactory.write_plot(hist_data,
                           f"{detector_name}_{scan_type}_dist_{l1}_{l2}",
                           out_format)


""" Plot the track position residual for the given variable """
def plot_track_pos_res(opts, detector, scan_type, plotFactory, out_format,
                       df1, label1, df2, label2, var):

    n_tracks = np.max(df1['track_id']) + 1
    tracks =  "rays" if scan_type == "ray" else "helices"
    # Where to place the legend box
    box_anchor_x = 1.02
    box_anchor_y = 1.28

    res = df1[var] - df2[var]

    # Remove outliers
    filter_res = np.absolute(res) < opts.outlier 
    filtered_res = res[filter_res]

    uOut = oOut = int(0)
    if not np.all(filter_res == True):
        print(f"\nRemoved outliers ({var}):")
        for i, r in enumerate(res):
            if math.fabs(r) > opts.outlier:
                track_id = (df1['track_id'].to_numpy())[i]
                print(f"track {track_id}: {df1[var][i]} - {df2[var][i]} = {r}")

                if r < 0.:
                    uOut = uOut + 1
                else:
                    oOut = oOut + 1


    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = plotting.legend_options('upper right', 4, 0.01, 0.0005)
    hist_data = plotFactory.hist1D(x        = filtered_res,
                                   figsize  = (9, 9),
                                   bins     = 100,
                                   xLabel   = r'$\mathrm{res}' + rf'\,{var}' + r'\,\mathrm{[mm]}$',
                                   setLog   = False,
                                   lgd_ops  = lgd_ops,
                                   uOutlier = uOut,
                                   oOutlier = oOut)

    mu, sig = plotFactory.fit_gaussian(hist_data)
    if mu is None or sig is None:
        print(rf"WARNING: fit failed (res ({tracks}): {label1} - {label2} )")

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    detector_name = detector.replace(' ', '_')
    l1 = label1.replace(' ', '_').replace("(", "").replace(")", "")
    l2 = label2.replace(' ', '_').replace("(", "").replace(")", "")

    plotFactory.write_plot(hist_data,
                           f"{detector_name}_{scan_type}_res_{var}_{l1}_{l2}",
                           out_format)
