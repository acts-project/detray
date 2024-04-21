# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# python includes
from collections import namedtuple
import math
import numpy as np

# python based plotting
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import matplotlib.style as style
style.use('tableau-colorblind10')
#style.use('seaborn-colorblind')

plt.rcParams.update({
    'text.usetex': True,
    'font.size': 20,
    'font.family': 'serif',
})

#-------------------------------------------------------------------------------
# Global identifiers
#-------------------------------------------------------------------------------

""" Pass plotting data between functions """
plt_data  = namedtuple('plt_data', 'fig ax lgd data bins mu rms errors')

""" Wrap the configuration for a legend """
legend_options = namedtuple('legend_options', 'loc ncol colspacing handletextpad')

""" Conveniently get the legend options """
def get_legend_options():
    return legend_options('upper right', 1, 1, 1)

#-------------------------------------------------------------------------------
# Data Plotting
#-------------------------------------------------------------------------------

"""
Plotter interface that uses pyplot/matplotlib.
"""
class pyplot_factory():

    def __init__(self, outDir, logger, atlas_badge = ""):
        self.name = 'Pyplot',
        self.outputPrefix = outDir
        self.logger = logger
        self.atlas_badge = atlas_badge
        self.badge_scale = 1.1


    """ Add legend to a plot. Labbels must be defined. """
    def add_legend(self, ax, options = get_legend_options()):
        return ax.legend(loc           = options.loc,
                         ncol          = options.ncol,
                         columnspacing = options.colspacing,
                         handletextpad = options.handletextpad)


    """
    Create a histogram from given input data. The normalization is achieved by 
    dividing the bin count by the total number of observations. The error is 
    calculated as the square root of the bin content.
    """
    def hist1D(self, x, w = None,
               xLabel = 'x', yLabel = '', title = "",  label = "",
               xMin   = None, xMax  = None, bins = 1,
               color  = 'tab:blue', alpha = 0.75,
               setLog    = False,
               normalize = False,
               showError = False,
               showStats = True,
               lgd_ops   = get_legend_options()):

        # Create fresh plot
        fig = plt.figure(figsize = (8, 6), layout='constrained')
        ax = fig.add_subplot(1, 1, 1)

        # Do calculations on data in the range of the histogram
        if not xMin is None and not xMax is None:
            x = x[np.where(x >= xMin)]
            x = x[np.where(x <= xMax)]
        else:
            xMin = np.min(x)
            xMax = np.max(x)

        # Nothing left to do
        if len(x) == 0:
            self.logger.debug(rf" create hist: empty data {label}")
            return plt_data(fig, ax, None, None, None, None, None, None)

        # Histogram normalization
        scale = 1./len(x) if normalize else 1.

        # Fill data
        data, bins, hist = ax.hist(x,
                                   weights   = w,
                                   range     = (xMin, xMax),
                                   bins      = bins,
                                   label     = f"{label}  ({len(x)} entries)",
                                   histtype  = 'stepfilled',
                                   density   = normalize,
                                   facecolor = mcolors.to_rgba(color, alpha),
                                   edgecolor = color)

        # Add some additional information
        if showStats:
            mean = np.mean(x, axis=0)
            rms  = np.sqrt(np.mean(np.square(x)))

            # Create empty plot with blank marker containing the extra label
            newline = '\n'
            ax.plot([], [], ' ', label= rf'mean = {mean:.2e}'
                                        rf'{newline}RMS  = {rms:.2e}')
        else:
            mean = None
            rms = None

        # Refine plot
        ax.set_title(title)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)
        ax.grid(True, alpha = 0.25)

        # Add legend
        lgd = self.add_legend(ax, lgd_ops)

        # Calculate the bin error
        binCenters = 0.5 * (bins[1:] + bins[:-1])
        errors     = np.sqrt(scale * data)
        if showError:
            ax.errorbar(binCenters, data,
                        yerr      = errors,
                        fmt       = '.',
                        linestyle = '',
                        color     = color)

        # Plot log scale
        if setLog:
            ax.set_yscale('log')

        return plt_data(fig, ax, lgd, data, bins, mean, rms, errors)


    """
    Create a 2D histogram from given input data. If z values are given they will
    be used as weights per bin.
    """
    def hist2D(self, x, y, z = None,
                xLabel = 'x', yLabel = 'y', zLabel = '', title = "", label = "",
                xMin   = None, xMax  = None, xBins = 1,
                yMin   = None, yMax  = None, yBins = 1,
                color  = 'tab:blue', alpha = 0.75,
                setLog    = False,
                showError = False,
                showStats = True):

        # Create fresh plot
        fig = plt.figure(figsize = (8, 6), layout='constrained')
        ax = fig.add_subplot(1, 1, 1)

        # Do calculations on data in the range of the histogram
        if not xMin is None and not xMax is None:
            x = x[np.where(x >= xMin)]
            x = x[np.where(x <= xMax)]
        else:
            xMin = np.min(x)
            xMax = np.max(x)

        if not yMin is None and not yMax is None:
            y = y[np.where(y >= yMin)]
            y = y[np.where(y <= yMax)]
        else:
            yMin = np.min(y)
            yMax = np.max(y)

        # Nothing left to do
        if len(x) == 0 or len(y) == 0:
            self.logger.debug(rf" create hist: empty data {label}")
            return plt_data(fig, ax, None, None, None, None, None, None)

        # Fill data
        data, xbins, ybins, hist = ax.hist2d(
                                    x, y, weights = z,
                                    range      = [(xMin, xMax), (yMin, yMax)],
                                    bins       = (xBins, yBins),
                                    label=f"{label}  ({len(x)*len(y)} entries)",
                                    facecolor  = mcolors.to_rgba(color, alpha),
                                    edgecolor  = None,
                                    rasterized = True)

        # Add some additional information
        if showStats:
            xMean = np.mean(x, axis=0)
            xRms  = np.sqrt(np.mean(np.square(x)))
            yMean = np.mean(y, axis=0)
            yRms  = np.sqrt(np.mean(np.square(y)))

            # Create empty plot with blank marker containing the extra label
            newline = '\n'
            ax.plot([], [], ' ', label= rf'xMean = {xMean:.2e}'
                                        rf'{newline}xRMS  = {xRms:.2e}'
                                        rf'yMean = {yMean:.2e}'
                                        rf'{newline}yRMS  = {yRms:.2e}')

        # Refine plot
        ax.set_title(title)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)

        # Add the colorbar
        fig.colorbar(hist, label = zLabel)

        return plt_data(fig, ax, None, data, None, None, None, None)


    """ Create a 2D scatter plot """
    def scatter(self, x, y,
                xLabel = "", yLabel = "", title = "", label = "",
                color  = 'tab:blue', alpha = 1,
                figsize   = (8, 6), 
                showStats = lambda x, _: f"{len(x)} entries", 
                lgd_ops   = get_legend_options()):

        fig = plt.figure(figsize = figsize, layout='constrained')
        ax = fig.add_subplot(1, 1, 1)

        # Create empty plot with blank marker containing the extra label
        ax.plot([], [], ' ', label=showStats(x, y))
        scatter = ax.scatter(x, y,
                             label = label,
                             c     = color,
                             s     = 0.1,
                             alpha = alpha,
                             rasterized=True)

        # Refine plot
        ax.set_title(title)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)
        ax.grid(True, alpha = 0.25)

        # Add legend
        lgd = self.add_legend(ax, lgd_ops)

        return plt_data(fig, ax, lgd, scatter, None, None, None, None)


    """ Add new data in a different color to a scatter plot """
    def highlight_region(self, plotData, x, y, color, label = ""):

        if label == "":
            plotData.ax.scatter(x, y, c = color, alpha = 1, s = 0.1,
                                rasterized=True)
        else:
            plotData.ax.scatter(x, y, c = color, alpha = 1, s = 0.1,
                                label=label, rasterized=True)

            # Update legend
            lgd = plotData.lgd
            handles, labels = lgd.axes.get_legend_handles_labels()
            lgd._legend_box = None
            lgd._init_legend_box(handles, labels)
            lgd._set_loc(lgd._loc)
            lgd.set_title(lgd.get_title().get_text())


    """ Safe a plot to disk """
    def write_plot(self, plot_data, name = "plot", file_format = "svg",
                   outPrefix = "", dpi = 450):
        if (outPrefix == ""):
            fileName = self.outputPrefix + name + "." + file_format
        else:
            fileName = outPrefix + name + "." + file_format

        plot_data.fig.savefig(fileName, dpi=dpi)
        plt.close(plot_data.fig)


    """ Safe a plot as svg """
    def write_svg(self, plot_data, name, outPrefix = ""):

        self.write_plot(plot_data, name, ".svg", outPrefix)


    """ Safe a plot as pdf """
    def write_pdf(self, plot_data, name, outPrefix = ""):

        self.write_plot(plot_data, name, ".pdf", outPrefix)


    """ Safe a plot as png """
    def write_png(self, plot_data, name, outPrefix = ""):

        self.write_plot(plot_data, name, ".png", outPrefix)
