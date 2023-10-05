# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# python includes
import numpy as np
import math
import matplotlib.pyplot as plt


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
def X0_vs_eta_phi(df, detector, plotFactory, format = "svg"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)

    # Plot the thickness of every material slab in units of X_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_tX0'],
                            title  = r'\textbf{Material Scan - thickness / }$\mathbf{X_0}$',
                            label  = detector,
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'$\mathbf{\phi}\,\mathrm{[rad]}$',
                            zLabel = r'thickness / $X_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "t_X0_map", format)

    # Plot path length through material of the respective ray in units of X_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_sX0'],
                            title  = r'\textbf{Material Scan - path length / }$\mathbf{X_0}$',
                            label  = detector,
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'$\mathbf{\phi}\,\mathrm{[rad]}$',
                            zLabel = r'path length / $X_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "s_X0_map", format)


""" Plot the material thickenss vs phi and eta in units of L_0 """
def L0_vs_eta_phi(df, detector, plotFactory, format = "svg"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)

    # Plot the thickness of every material slab in units of L_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_tL0'],
                            title  = r'\textbf{Material Scan - thickness /} $\mathbf{\Lambda_0}$',
                            label  = detector,
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'$\mathbf{\phi}\,\mathrm{[rad]}$',
                            zLabel = r'thickness / $\Lambda_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "t_L0_map", format)

    # Plot path length through material of the respective ray in units of L_0
    hist_data = plotFactory.hist2D(
                            x      = df["eta"],
                            y      = df["phi"],
                            z      = df['mat_sL0'],
                            title  = r'\textbf{Material Scan - path length /} $\mathbf{\Lambda_0}$',
                            label  = detector,
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'$\mathbf{\phi}\,\mathrm{[rad]}$',
                            zLabel = r'path length / $\Lambda_0$',
                            xBins = xBinning, yBins = yBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "s_L0_map", format)


""" Plot the material thickness in units of X_0 vs eta """
def X0_vs_eta(df, detector, plotFactory, format = "svg"):

    # Histogram bin edges
    xBinning, _ = get_n_bins(df)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_tX0'],
                            normalize = False,
                            title  = r'\textbf{Material Scan - thickness / }$\mathbf{X_0}$',
                            label  = rf'{detector}',
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'thickness / $X_0$',
                            bins = xBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "t_X0", format)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_sX0'],
                            normalize = False,
                            title  = r'\textbf{Material Scan - path length / }$\mathbf{X_0}$',
                            label  = rf'{detector}',
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'path length / $X_0$',
                            bins = xBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "s_X0", format)


""" Plot the material thickness in units of L_0 vs eta """
def L0_vs_eta(df, detector, plotFactory, format = "svg"):

    # Histogram bin edges
    xBinning, _ = get_n_bins(df)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_tL0'],
                            normalize = False,
                            title  = r'\textbf{Material Scan - thickness / }$\mathbf{\Lambda_0}$',
                            label  = rf'{detector}',
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'thickness / $\Lambda_0$',
                            bins = xBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "t_L0", format)

    hist_data = plotFactory.hist1D(
                            x      = df['eta'],
                            w      = df['mat_sL0'],
                            normalize = False,
                            title  = r'\textbf{Material Scan - path length / }$\mathbf{\Lambda_0}$',
                            label  = rf'{detector}',
                            xLabel = r'$\mathbf{\eta}$',
                            yLabel = r'path length / $\Lambda_0$',
                            bins = xBinning,
                            showStats = False)

    plotFactory.write_plot(hist_data, "s_L0", format)
