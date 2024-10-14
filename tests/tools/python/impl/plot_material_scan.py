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
def read_material_data(inputdir, logging, det_name, read_cuda):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    detector_name = "default_detector"
    material_scan_file = cpu_material_trace_file = cuda_material_trace_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if filename.find(det_name + '_material_scan') != -1:
            material_scan_file = inputdir + "/" + filename
        elif filename.find(det_name + '_cpu_navigation_material_trace') != -1:
            cpu_material_trace_file = inputdir + "/" + filename
        elif read_cuda and filename.find(det_name + '_cuda_navigation_material_trace') != -1:
            cuda_material_trace_file = inputdir + "/" + filename

    df_scan = pd.read_csv(material_scan_file,
                               float_precision='round_trip')
    df_cpu_trace = pd.read_csv(cpu_material_trace_file,
                               float_precision='round_trip')
    df_cuda_trace = pd.DataFrame({})
    if read_cuda:
        df_cuda_trace = pd.read_csv(cuda_material_trace_file,
                               float_precision='round_trip')

    return df_scan, df_cpu_trace, df_cuda_trace


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


""" Calculate the binwise errors: Standard Error on the Mean """
def get_errors(df, n, name):
    # Number of entries per bin
    errors = []
    # Iterate over data per bin
    for i in range(0, n):
        # Project the next n rows of the data frame (the bin content is the 
        # mean of the material along phi)
        bin_data = df.iloc[i * n:(i + 1) * n]
        # Caluculate the error on the mean: stddev/sqrt(n)
        errors.append(np.std(bin_data[name], axis=0) / math.sqrt(n))

    return errors


""" Plot the material thickenss vs phi and eta in units of X_0 """
def X0_vs_eta_phi(df, label, detector, plotFactory,  out_format =  "pdf"):

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
                            figsize  = (9, 7),
                            showStats = False)

    plotFactory.write_plot(hist_data, detector + "_" + label + "_t_X0_map",  out_format)

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
                            figsize  = (9, 7),
                            showStats = False)

    plotFactory.write_plot(hist_data, detector + "_" + label + "_s_X0_map",  out_format)


""" Plot the material thickenss vs phi and eta in units of L_0 """
def L0_vs_eta_phi(df, label, detector, plotFactory,  out_format =  "pdf"):

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
                            figsize  = (9, 7),
                            showStats = False)

    plotFactory.write_plot(hist_data, detector + "_" + label + "_t_L0_map",  out_format)

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
                            figsize  = (9, 7),
                            showStats = False)

    plotFactory.write_plot(hist_data, detector + "_" + label + "_s_L0_map",  out_format)


""" Plot the material thickness in units of X_0 vs eta """
def X0_vs_eta(df, label, detector, plotFactory,  out_format =  "pdf"):
    # Where to place the legend box
    box_anchor_x = 1.02
    box_anchor_y = 1.145

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)
    lgd_ops = plotting.get_legend_options()
    lgd_ops._replace(loc = 'upper center')

    # Same number of entries in every bin as per uniform ray scan
    n_phi = len(yBinning) - 1

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_tX0'] / n_phi,
                            errors = get_errors(df, n_phi, 'mat_tX0'),
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'thickness / $X_0$',
                            bins = xBinning,
                            showStats = False,
                            figsize  = (9, 7),
                            lgd_ops = lgd_ops)

    # Move the legend ouside plot
    hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(hist_data, detector + "_" + label + "_t_X0",  out_format)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_sX0'] / n_phi,
                            errors = get_errors(df, n_phi, 'mat_sX0'),
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'path length / $X_0$',
                            bins = xBinning,
                            showStats = False,
                            figsize  = (9, 7),
                            lgd_ops = lgd_ops)

    # Move the legend ouside plot
    hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(hist_data, detector + "_" + label + "_s_X0",  out_format)


""" Plot the material thickness in units of L_0 vs eta """
def L0_vs_eta(df, label, detector, plotFactory,  out_format =  "pdf"):
    # Where to place the legend box
    box_anchor_x = 1.02
    box_anchor_y = 1.145

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)
    lgd_ops = plotting.get_legend_options()
    lgd_ops._replace(loc = 'upper center')

    # Same number of entries in every bin as per uniform ray scan
    n_phi = len(yBinning) - 1

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_tL0'] / n_phi,
                            errors = get_errors(df, n_phi, 'mat_tL0'),
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'thickness / $\Lambda_0$',
                            bins = xBinning,
                            showStats = False,
                            figsize  = (9, 7),
                            lgd_ops = lgd_ops)

    # Move the legend ouside plot
    hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(hist_data, detector + "_" + label + "_t_L0",  out_format)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_sL0'] / n_phi,
                            errors = get_errors(df, n_phi, 'mat_sL0'),
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = r'$\eta$',
                            yLabel = r'path length / $\Lambda_0$',
                            bins = xBinning,
                            showStats = False,
                            figsize  = (9, 7),
                            lgd_ops = lgd_ops)

    # Move the legend ouside plot
    hist_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    plotFactory.write_plot(hist_data, detector + "_" + label + "_s_L0",  out_format)


""" Compare two material distributions """
def compare_mat(df_truth, df_rec, label, detector, plotFactory, out_format =  "pdf"):
    # Where to place the legend box
    box_anchor_x = 1.02
    box_anchor_y = 1.29

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df_truth)
    lgd_ops = plotting.get_legend_options()
    lgd_ops._replace(loc = 'upper center')

    # Same number of entries in every bin as per uniform ray scan
    n_phi = len(yBinning) - 1

    truth_data = plotFactory.hist1D(
                            x      = df_truth['eta'],
                            w      = df_truth['mat_sX0'] / n_phi,
                            errors = get_errors(df_truth, n_phi, 'mat_sX0'),
                            normalize = False,
                            label  = rf'{detector}: scan',
                            xLabel = r'$\eta$',
                            yLabel = r'path length / $X_0$',
                            bins = xBinning,
                            showStats = False,
                            figsize  = (10, 10),
                            layout = 'tight',
                            lgd_ops = lgd_ops)

    # Move the legend ouside plot
    truth_data.lgd.set_bbox_to_anchor((box_anchor_x, box_anchor_y))

    # Add recorded data for comparison
    rec_data = plotFactory.add_plot(
                            oldHist   = truth_data,
                            x         = df_rec['eta'].to_numpy(dtype = np.double),
                            w         = df_rec['mat_sX0'] / n_phi,
                            errors    = get_errors(df_rec, n_phi, 'mat_sX0'),
                            normalize = False,
                            label     = rf'{detector}: navigator')

    # Add a ratio plot to hist_data
    ratio_data = plotFactory.add_ratio(
                        nom        = truth_data,
                        denom      = rec_data,
                        label      = 'scan/navigation',
                        setLog     = False,
                        showErrors = True)

    plotFactory.write_plot(ratio_data, detector + "_" + label + "_comparison_s_X0",  out_format)
