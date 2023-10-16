# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import plot_helpers
from pyplot_factory import legend_options

# python includes
import numpy as np

""" Plot the intersection points of the detector with the rays - xy view """
def intersection_points_xy(opts, df, detector, scan_type, plotFactory,  out_format = "png"):

    n_rays = np.max(df['index']) + 1
    tracks = "rays" if scan_type == "ray" else "helices"

    # Reduce data to the requested z-range (50mm tolerance)
    min_z = opts.z_range[0]
    max_z = opts.z_range[1]
    assert min_z < max_z, "xy plotting range: min z must be smaller that max z"
    sensitive_range = lambda data: ((data['z'] > min_z) & (data['z'] < max_z) & (data['type'] == 1))
    portal_range = lambda data: ((data['z'] > min_z) & (data['z'] < max_z) & (data['type'] == 0))
    passive_range = lambda data: ((data['z'] > min_z) & (data['z'] < max_z) & (data['type'] == 2))

    senstive_x, senstive_y = plot_helpers.filter_data(
                                                data      = df,
                                                filter    = sensitive_range,
                                                variables = ['x', 'y'])

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = legend_options('upper center', 4, 0.4, 0.005)
    hist_data = plotFactory.scatter(
                            figsize = (8, 8),
                            x       = senstive_x,
                            y       = senstive_y,
                            xLabel  = r'$x\,\mathrm{[mm]}$',
                            yLabel  = r'$y\,\mathrm{[mm]}$',
                            label   = "sensitives",
                            color   = 'C5',
                            showStats = lambda x, y: f"{n_rays} {tracks}",
                            lgd_ops = lgd_ops)

    # Portal surfaces
    if not opts.hide_portals:
        portal_x, portal_y = plot_helpers.filter_data(
                                                data      = df,
                                                filter    = portal_range,
                                                variables = ['x', 'y'])

        plotFactory.highlight_region(hist_data, portal_x, portal_y, 'C0',      \
                                     "portals")

    # Passive surfaces
    if not opts.hide_passives:
        passive_x, passive_y = plot_helpers.filter_data(
                                                data      = df,
                                                filter    = passive_range,
                                                variables = ['x', 'y'])

        plotFactory.highlight_region(hist_data, passive_x, passive_y, 'C2',    \
                                     "passives")

    # Refine legend
    hist_data.lgd.legendHandles[0].set_visible(False)
    for handle in hist_data.lgd.legendHandles[1:]:
        handle.set_sizes([40])

    # For this plot, move the legend ouside
    hist_data.lgd.set_bbox_to_anchor((0.5, 1.11))

    # Adjust spacing in box
    for vpack in hist_data.lgd._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

    detector_name = detector.replace(' ', '_')
    plotFactory.write_plot(hist_data, f"{detector_name}_{scan_type}_scan_xy",  out_format)


""" Plot the intersection points of the detector with the rays - rz view """
def intersection_points_rz(opts, df, detector, scan_type, plotFactory,  out_format = "png"):

    n_rays = np.max(df['index']) + 1
    tracks =  "rays" if scan_type == "ray" else "helices"

    # Reduce data to the requested z-range
    sensitive_range = lambda data: (data['type'] == 1)
    portal_range = lambda data: (data['type'] == 0)
    passive_range = lambda data: (data['type'] == 2)

    sensitive_x, sensitive_y, sensitive_z = plot_helpers.filter_data(
                                                data      = df,
                                                filter    = sensitive_range,
                                                variables = ['x', 'y', 'z'])

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = legend_options('upper center', 4, 0.8, 0.005)
    hist_data = plotFactory.scatter(
                            figsize = (12, 6),
                            x      = sensitive_z,
                            y      = np.hypot(sensitive_x, sensitive_y),
                            xLabel = r'$z\,\mathrm{[mm]}$',
                            yLabel = r'$r\,\mathrm{[mm]}$',
                            label  = "sensitives",
                            color  = 'C5',
                            showStats = lambda x, y: f"{n_rays} {tracks}",
                            lgd_ops = lgd_ops)

    # Portal surfaces
    if not opts.hide_portals:
        portal_x, portal_y, portal_z = plot_helpers.filter_data(
                                                    data      = df,
                                                    filter    = portal_range,
                                                    variables = ['x', 'y', 'z'])

        plotFactory.highlight_region(hist_data, portal_z,                      \
                                     np.hypot(portal_x, portal_y),             \
                                     'C0', "portals")

    # Passive surfaces
    if not opts.hide_passives:
        passive_x, passive_y, passive_z = plot_helpers.filter_data(
                                                    data      = df,
                                                    filter    = passive_range,
                                                    variables = ['x', 'y', 'z'])

        plotFactory.highlight_region(hist_data, passive_z,                     \
                                    np.hypot(passive_x, passive_y),            \
                                    'C2', "passives")

    # Refine legend
    hist_data.lgd.legendHandles[0].set_visible(False)
    for handle in hist_data.lgd.legendHandles[1:]:
        handle.set_sizes([40])

    # For this plot, move the legend ouside
    hist_data.lgd.set_bbox_to_anchor((0.5, 1.15))

    # Adjust spacing in box
    for vpack in hist_data.lgd._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

    detector_name = detector.replace(' ', '_')
    plotFactory.write_plot(hist_data, f"{detector_name}_{scan_type}_scan_rz",  out_format)
