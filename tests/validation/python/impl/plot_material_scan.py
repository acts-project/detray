# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
import plotting

# python includes
import math
import numpy as np
import os
import pandas as pd


""" Read the material scan data from file and prepare data frame """
def read_material_data(inputdir, logging):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    detector_name = "default_detector"
    material_scan_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if filename.find('material_scan_') != -1:
            material_scan_file = inputdir + "/" + filename
            file_name = os.path.basename(material_scan_file)
            detector_name = file_name.removeprefix('material_scan_ ')
            detector_name = detector_name.removesuffix('.csv')

    detector_name = detector_name.replace('_', ' ')

    df = pd.read_csv(material_scan_file)

    return detector_name, df


""" Calculate edges of bins to plot the mateiral data """
def get_n_bins(df):
    # Find the number of ray directions
    row_count = df.groupby(df['eta']).count()
    yBins = row_count['phi'].max()
    xBins = int(len(df['eta']) / yBins)
    assert len(df['eta']) == xBins * yBins, "Could not infer the number of rays correctly"

    # Get the axis spacing
    x_range = np.max(df['eta']) - np.min(df['eta'])
    xBinning = np.linspace(np.min(df['eta']) - 0.5 * x_range/xBins, np.max(df['eta']) + 0.5 * x_range/xBins, xBins + 1)

    y_range = np.max(df['phi']) - np.min(df['phi'])
    yBinning = np.linspace(np.min(df['phi']) - 0.5 * y_range/yBins, np.max(df['phi']) + 0.5 * y_range/yBins, yBins + 1)

    return xBinning, yBinning


""" Plot the material thickenss vs phi and eta in units of X_0 """
def X0_vs_eta_phi(df, detector, plotFactory,  out_format =  "pdf"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)

    # Plot the thickness of every material slab in units of X_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_tX0'],
                            label  = detector,
                            xLabel = r'$\eta$',
                            yLabel = r'$\phi\,\mathrm{[rad]}$',
                            zLabel = r'thickness / $X_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "t_X0_map",  out_format)

    # Plot path length through material of the respective ray in units of X_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_sX0'],
                            label  = detector,
                            xLabel = r'$\eta$',
                            yLabel = r'$\phi\,\mathrm{[rad]}$',
                            zLabel = r'path length / $X_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "s_X0_map",  out_format)


""" Plot the material thickenss vs phi and eta in units of L_0 """
def L0_vs_eta_phi(df, detector, plotFactory,  out_format =  "pdf"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)

    # Plot the thickness of every material slab in units of L_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_tL0'],
                            label  = detector,
                            xLabel = r'$\eta$',
                            yLabel = r'$\phi\,\mathrm{[rad]}$',
                            zLabel = r'thickness / $\Lambda_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "t_L0_map",  out_format)

    # Plot path length through material of the respective ray in units of L_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_sL0'],
                            label  = detector,
                            xLabel = r'$\eta$',
                            yLabel = r'$\phi\,\mathrm{[rad]}$',
                            zLabel = r'path length / $\Lambda_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "s_L0_map",  out_format)


""" Plot the material thickness in units of X_0 vs eta """
def X0_vs_eta(df, detector, plotFactory,  out_format =  "pdf"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)
    lgd_ops = plotting.get_legend_options()
    lgd_ops._replace(loc = 'upper center')

    phi_normalization = len(yBinning) - 1

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_tX0'] / phi_normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'thickness / $X_0$',
                            bins = xBinning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, "t_X0",  out_format)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_sX0'] / phi_normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'path length / $X_0$',
                            bins = xBinning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, "s_X0",  out_format)


""" Plot the material thickness in units of L_0 vs eta """
def L0_vs_eta(df, detector, plotFactory,  out_format =  "pdf"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)
    lgd_ops = plotting.get_legend_options()
    lgd_ops._replace(loc = 'upper center')

    phi_normalization = len(yBinning) - 1

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_tL0'] / phi_normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'thickness / $\Lambda_0$',
                            bins = xBinning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, "t_L0",  out_format)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_sL0'] / phi_normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'path length / $\Lambda_0$',
                            bins = xBinning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, "s_L0",  out_format)
