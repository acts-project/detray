# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
import plotting

# python includes
from collections import namedtuple
import json
import itertools
import os
import pandas as pd
import sys

# How to label plots
label_data = namedtuple("label_data", "title label x_axis y_axis")

# Define labels for google benchmark data collections
label_dict = {
    "real_time": label_data(
        "Propagation Latency", "", "No. tracks", r"t $[\mathrm{ms}]$"
    ),
    "TracksPropagated": label_data(
        "Propagation Throughout", "", "No. tracks", r"Prop. rate $[\mathrm{MHz}]$"
    ),
}

# Common options
ldg_loc = "upper left"


""" Read google benchmark data from json file """


def read_benchmark_data(logging, input_path, benchmark_file):

    file_path = input_path + benchmark_file
    with open(file_path, "r") as file:
        logging.debug(f"Reading file '{file_path}'")

        results = json.load(file)

        context = results["context"]
        data = pd.DataFrame(results["benchmarks"])

        return context, data

    logging.error(f"Could not find file: {benchmark_file}")

    return None, None


""" Adds a column 'x' to the data frame that contains the number of tracks """


def add_track_multiplicity_column(df):

    assert "_TRACKS" in str(
        df["run_name"][0]
    ), "Benchmark case name not correctly formatted: (BM_PROPAGATION_<detector name>_<#tracks>_TRACKS)"

    # The number of tracks is the second last part of the benchmark name
    find_track_multiplicity = lambda n: (int(n.split("_")[-2]))

    # Add new column based on benchmark case name
    df["x"] = df["run_name"].apply(find_track_multiplicity)


""" Read the benchmark data and prepare it for plotting """


def prepare_benchmark_data(logging, input_dir, file):

    # Convert benchmark timings to 'ms'
    unit_conversion = {"ns": 10**-6, "um": 10**-3, "ms": 1, "s": 10**3}

    # Read the data part into a pandas frame
    context, data = read_benchmark_data(logging, input_dir, file)

    if context is None or data is None:
        logging.warning(f"Failed to read data in file: {file}")
        sys.exit(1)

    # Add the number of tracks per benchmark case as new column 'x'
    # A column called 'x' is expected by the 'plot_benchmark' method
    add_track_multiplicity_column(data)

    # Convert timings to 'ms'
    bench_time_unit = data["time_unit"][0]
    to_milliseconds = lambda x: (x * unit_conversion[bench_time_unit])

    data["real_time"] = data["real_time"].apply(to_milliseconds)
    data["cpu_time"] = data["cpu_time"].apply(to_milliseconds)

    # Convert from Hz to MHz
    data["TracksPropagated"] = data["TracksPropagated"] / 1000000

    return context, data


"""
Plot the benchmark latency and throughout for different hardware backends and
algebra plugins
"""


def plot_benchmark_case(
    context,
    df,
    plot_factory,
    label,
    title="",
    data_type="real_time",
    marker=".",
    plot=None,
    yaxis_format=None,
):

    assert len(df["x"]) != 0, "Data frame has to provide column 'x'"
    assert len(df[data_type]) != 0, f"Data frame has to provide column '{data_type}'"

    # Filter the relevant data from the frame
    median = lambda data_frame: (data_frame["aggregate_name"] == "median")
    stddev = lambda data_frame: (data_frame["aggregate_name"] == "stddev")

    data, n_tracks = plotting.filter_data(
        data=df, filter=median, variables=[data_type, "x"]
    )

    stddev = plotting.filter_data(data=df, filter=stddev, variables=[data_type])

    if plot is None:
        # Create new plot
        lgd_ops = plotting.legend_options(
            loc=ldg_loc, horiz_anchor=1.0, vert_anchor=1.02
        )

        labels = label_dict[data_type]
        x_axis_opts = plotting.axis_options(
            label=labels.x_axis, log_scale=True, tick_positions=n_tracks
        )
        y_axis_opts = plotting.axis_options(
            label=labels.y_axis, log_scale=True, label_format=yaxis_format
        )

        # Plot the propagation latency against the number of tracks
        plot_data = plot_factory.graph(
            x=n_tracks,
            y=data,
            y_errors=stddev,
            x_axis=x_axis_opts,
            y_axis=y_axis_opts,
            title=title,
            label=label,
            lgd_ops=lgd_ops,
            marker=marker,
            figsize=(15, 8),
        )
    else:
        # Add new data to exiting plot
        plot_data = plot_factory.add_graph(
            plot=plot,
            x=n_tracks,
            y=data,
            y_errors=stddev,
            label=label,
            marker=marker,
            color=None,
        )

    return plot_data


""" Plot the data of all benchmark files given in 'data_files' """


def plot_benchmark_data(
    logging,
    input_dir,
    det_name,
    file_list,
    label_list,
    title,
    plot_series_name,
    plot_factory,
    out_format,
):

    # Cylce through marker styles per plot
    marker_styles = ["o", "x", "*", "v", "s", "^", "<", ">"]

    # Plot types for benchmarks
    benchmark_plots = namedtuple("benchmark_plots", "latency throughput")

    # Save the different plots per hardware backend
    plots = benchmark_plots(None, None)
    marker_style_cycle = itertools.cycle(marker_styles)

    # Go through all benchmark data files in the list and make a comparison plot
    for i, file in enumerate(file_list):
        # Get the data for the next benchmark case
        context, data = prepare_benchmark_data(logging, input_dir, file)
        marker = next(marker_style_cycle)

        # Initialize plots
        if i == 0:

            # Plot the data against the number of tracks
            latency_plot = plot_benchmark_case(
                context=context,
                df=data,
                plot_factory=plot_factory,
                label=label_list[i],
                data_type="real_time",
                marker=marker,
                title=title,
                yaxis_format=None,
            )

            throughput_plot = plot_benchmark_case(
                context=context,
                df=data,
                plot_factory=plot_factory,
                label=label_list[i],
                data_type="TracksPropagated",
                marker=marker,
                title=title,
                yaxis_format="{x:3.2f}",
            )

            plots = benchmark_plots(latency_plot, throughput_plot)

        # Add new data to plots
        else:
            plot_benchmark_case(
                context=context,
                df=data,
                plot_factory=plot_factory,
                label=label_list[i],
                data_type="real_time",
                marker=marker,
                plot=plots.latency,
            )

            plot_benchmark_case(
                context=context,
                df=data,
                plot_factory=plot_factory,
                label=label_list[i],
                data_type="TracksPropagated",
                marker=marker,
                plot=plots.throughput,
            )

    # Write to disk
    plot_factory.write_plot(
        plots.latency, f"{det_name}_{plot_series_name}_latency", out_format
    )

    plot_factory.write_plot(
        plots.throughput, f"{det_name}_{plot_series_name}_throughput", out_format
    )
