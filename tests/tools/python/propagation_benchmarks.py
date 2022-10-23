# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray imports
from impl import read_benchmark_data, plot_benchmark_data
from options import (
    common_options,
    detector_io_options,
    random_track_generator_options,
    propagation_options,
    plotting_options,
)
from options import (
    parse_common_options,
    parse_detector_io_options,
    parse_plotting_options,
)
from plotting import pyplot_factory as plt_factory
from utils import read_detector_name
from utils import add_track_generator_args, add_propagation_args, add_detector_io_args

# python imports
import argparse
from collections import namedtuple
import os
import platform
import subprocess
import sys

# Known hardware backend types
bknd_types = ["cpu", "cuda", "sycl"]

# Patterns to be removed from processor names for simplicity
bknd_patterns = ["CPU", "(TM)", "GHz", "@"]


# Simpler hardware backend tag
def __compactify_bknd_name(name, patterns=bknd_patterns):
    out = ""
    for sub_string in name.split(" "):
        if any(p in sub_string for p in patterns):
            continue

        out = f"{out} {sub_string}"

    # Remove preceeding whitespace
    return out[1:]


# Peek into the benchmark context to get the name of the backend
def __read_backend_name(logging, input_dir, bknd, data_file):
    context, _ = read_benchmark_data(logging, input_dir, data_file)
    bknd_name = __compactify_bknd_name(context["CPU" if "cpu" in bknd else "GPU"])

    return bknd_name


# Parse and check the user provided input data files
def __parse_input_data_files(args, det_name, algebra_plugins):
    input_data_files = []
    for file in args.data_files:
        if not os.path.isfile(file):
            logging.error(f"File not found! ({file})")
            sys.exit(1)

        file_name, file_extension = os.path.splitext(file)
        format_msg = f"Benchmark data file name needs to be of the form <detector>_benchmark_data_<cpu|cuda|sycl>_<algebra-plugin>.json: e.g. 'toy_detector_benchmark_data_cpu_eigen.json' ({file})"

        if f"{det_name}_benchmark_data" not in file_name:
            logging.error("Wrong prefix: " + format_msg)
            sys.exit(1)
        if file_extension != ".json":
            logging.error("Wrong file extension. Should be '.json': " + format_msg)
            sys.exit(1)
        if not any(p in file_name for p in bknd_types):
            logging.error(
                "No hardware backend type found (cpu|cuda|sycl): " + format_msg
            )
            sys.exit(1)
        if not any(p in file_name for p in algebra_plugins):
            logging.error("No algebra-plugin name found: " + format_msg)
            sys.exit(1)

        input_data_files.append(file)

    return input_data_files


# Gather and check benchmark executables and resulting data files for every
# hardware backend type and algebra plugin
def __generate_benchmark_dict(
    args, logging, bindir, det_name, input_data_files, algebra_plugins
):
    benchmark_files = namedtuple("benchmark_files", "bin data_files")
    benchmarks = {"cpu": benchmark_files([], [])}
    if args.cuda:
        benchmarks["cuda"] = benchmark_files([], [])
    if args.sycl:
        # benchmarks["sycl"] = benchmark_files([], [])
        logging.error("SYCL propagation benchmark is not implemented")

    for bknd, files in benchmarks.items():
        for algebra in algebra_plugins:
            binary = f"{bindir}/detray_propagation_benchmark_{bknd}_{algebra}"
            data_file = f"{det_name}_benchmark_data_{bknd}_{algebra}.json"

            # If the results should not be read from file, run the benchmark
            if data_file not in (os.path.basename(f) for f in input_data_files):
                # Register binary if it exists
                if not os.path.isdir(bindir) or not os.path.isfile(binary):
                    logging.warning(
                        f"Propagation benchmark binary not found! ({binary})"
                    )
                else:
                    files.bin.append(binary)
                    files.data_files.append(data_file)
            else:
                for f in input_data_files:
                    if data_file == os.path.basename(f):
                        # Add result file with custom path to be plotted
                        files.data_files.append(f)
    return benchmarks


def __main__():

    # ---------------------------------------------------------------arg parsing

    descr = "Detray Propagation Benchmark"

    # Define options
    parent_parsers = [
        common_options(descr),
        detector_io_options(),
        random_track_generator_options(),
        propagation_options(),
        plotting_options(),
    ]

    parser = argparse.ArgumentParser(description=descr, parents=parent_parsers)

    parser.add_argument(
        "--bindir",
        "-bin",
        help=("Directoy containing the benchmark executables"),
        default="./bin",
        type=str,
    )
    parser.add_argument(
        "--cuda",
        help=("Run the CUDA propagation benchmarks."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--sycl",
        help=("Run the SYCL propagation benchmarks (Not implemented)."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--sort_tracks",
        help=("Sort the track samples by theta."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--interleave",
        help=("Interleave the benchmark cases randomly."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--benchmark_repetitions",
        help=("Number of repeated benchmark runs."),
        default=2,
        type=int,
    )
    parser.add_argument(
        "--algebra_plugins",
        "-ap",
        nargs="*",
        help=(
            "Algebra plugins to be benchmarked (the plugin must be enabled at build time)."
        ),
        default=["array"],
        type=str,
    )
    parser.add_argument(
        "--data_files",
        "-f",
        nargs="*",
        help=("Read the benchmark results from a Google benchmark json file instead."),
        default=[],
        type=str,
    )

    # Parse options
    args = parser.parse_args()

    logging = parse_common_options(args, descr)
    parse_detector_io_options(args, logging)
    input_dir, out_dir, out_format = parse_plotting_options(args, logging)

    # Check bin path
    bindir = args.bindir.strip("/")

    # Get detector name
    det_name = read_detector_name(args.geometry_file, logging)
    logging.debug("Detector: " + det_name)

    # Unique set of algebra plugins to be included in the plots
    algebra_plugins = set(args.algebra_plugins)

    # Check user provided benchmark result files
    input_data_files = __parse_input_data_files(args, det_name, algebra_plugins)

    # Get dictionary of benchmark files per hardware backend type
    benchmarks = __generate_benchmark_dict(
        args, logging, bindir, det_name, input_data_files, algebra_plugins
    )

    # -----------------------------------------------------------------------run

    # Pass on the options for the detray benchmark executable
    args_list = []

    # Add parsed options to argument list
    add_detector_io_args(args_list, args)
    add_track_generator_args(args_list, args)
    add_propagation_args(args_list, args)

    if args.sort_tracks:
        args_list.append("--sort_tracks")

    logging.debug(args_list)

    # Pass on the options for google benchmark
    benchmark_options = [
        f"--benchmark_repetitions={args.benchmark_repetitions}",
        # "--benchmark_min_time=50x", taken from user guide, but does not work...
        "--benchmark_display_aggregates_only=true",
        # "--benchmark_time_unit=ms", taken from user guide, but does not work...
        "--benchmark_out_format=json",
        f"--benchmark_context=CPU={platform.processor()}",
    ]
    if args.interleave:
        benchmark_options.append("--benchmark_enable_random_interleaving=true")

    # Run the benchmarks
    for bknd, files in benchmarks.items():

        if args.cuda or args.sycl:
            # Try to get the GPU name
            gpu = ""
            try:
                gpu = str(subprocess.check_output(["nvidia-smi", "-L"]))
            except:
                gpu = "Unknown"

            benchmark_options.append(f"--benchmark_context=GPU={gpu}")

        for binary in files.bin:
            algebra = binary.split(f"benchmark_{bknd}_")[-1]
            subprocess.run(
                [
                    binary,
                    f"--benchmark_context=Algebra={algebra}",
                    f"--benchmark_out=./{det_name}_benchmark_data_{bknd}_{algebra}.json",
                ]
                + benchmark_options
                + args_list
            )

    # ----------------------------------------------------------------------plot

    logging.info("Generating plots...\n")

    plot_factory = plt_factory(out_dir, logging)

    # Plot all data files per hardware backend
    # (comparison of different algebra-plugins)
    for bknd, benchmark_data in benchmarks.items():
        bknd_name = __read_backend_name(
            logging, input_dir, bknd, benchmark_data.data_files[0]
        )

        # Generate plot labels
        plot_labels = []
        for file in benchmark_data.data_files:
            # Get algebra-plugin from file name
            file_stem, _ = os.path.splitext(file)
            algebra = file_stem.split(f"{det_name}_benchmark_data_{bknd}_")[-1]
            plot_labels.append(algebra)

        plot_benchmark_data(
            logging,
            input_dir,
            det_name,
            benchmark_data.data_files,
            plot_labels,
            f"hardware backend: {bknd.upper()} ({bknd_name})",
            f"prop_benchmark_algebra-plugin_comparison_{bknd}",
            plot_factory,
            out_format,
        )

    # Plot results for different hardware backends using the same algebra plugin
    # (comparison of different hardware backends)
    for algebra in algebra_plugins:
        data_files_per_plugin = []
        plot_labels = []

        for bknd, benchmark_data in benchmarks.items():
            bknd_name = __read_backend_name(
                logging, input_dir, bknd, benchmark_data.data_files[0]
            )

            for data_file in benchmark_data.data_files:
                if algebra in data_file:
                    data_files_per_plugin.append(data_file)
                    plot_labels.append(f"{bknd}: {bknd_name}")

        plot_benchmark_data(
            logging,
            input_dir,
            det_name,
            data_files_per_plugin,
            plot_labels,
            f"algebra-plugin: {algebra}",
            f"prop_benchmark_backend_comparison_{algebra}",
            plot_factory,
            out_format,
        )


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
