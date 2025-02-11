# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray imports
from impl import read_benchmark_data, plot_benchmark_data, plot_scaling_data
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
bknd_patterns = [
    "CPU",
    "(TM)",
    "GHz",
    "@",
    "Core",
    "Processor",
    "with",
    "Radeon",
    "Graphics",
    "GB",
    "Laptop",
    "GPU",
    "GeForce",
]


# Simpler hardware backend tag
def __compactify_bknd_name(name, patterns=bknd_patterns):
    out = ""
    for sub_string in name.split(" "):
        if any(p in sub_string for p in patterns):
            continue

        out = f"{out} {sub_string}"

    # Remove preceeding whitespace
    return name if len(out[1:]) == 0 else out[1:]


# Peek into the benchmark context to get the name of the backend
def __read_context_metadata(logging, input_dir, data_file):
    context, _ = read_benchmark_data(logging, input_dir, data_file)
    bknd = context["Backend"]
    bknd_name = __compactify_bknd_name(context["Backend Name"])
    algebra = context["Algebra-plugin"]
    setup = context["Detector Setup"]

    return bknd, bknd_name, algebra, setup


# Parse and check the user provided input data files
def __parse_input_data_files(args, logging):
    input_data_files = []
    for file in args.data_files:
        if not os.path.isfile(file):
            logging.error(f"File not found! ({file})")
            sys.exit(1)

        _, file_extension = os.path.splitext(file)

        if file_extension != ".json":
            logging.error("Wrong file extension. Should be '.json': " + format_msg)
            sys.exit(1)

        input_data_files.append(file)

    return input_data_files


# Gather and check benchmark executables and resulting data files for every
# hardware backend type and algebra plugin
def __generate_benchmark_dict(
    args,
    logging,
    bindir,
    det_name,
    input_data_files,
    algebra_plugins,
    bench_type="benchmark",
):
    # Bundle benchmark metadata
    benchmark_metadata = namedtuple(
        "benchmark_metadata",
        "name algebra setup bin file",
        defaults=["Unknown", "Unknown Algebra", "", None, None],
    )

    # Resulting dictionary
    benchmarks = {"CPU": {}}
    if args.cuda:
        benchmarks["CUDA"] = {}
    if args.sycl:
        # benchmarks["SYCL"] = {}
        logging.error(f"SYCL propagation {bench_type} is not implemented")

    # Register the input data files for plotting
    for f in input_data_files:
        if not os.path.isfile(f):
            logging.error(f"File does not exist: {f}")
            sys.exit(1)

        # Add file to benchmark dict
        input_dir = os.path.dirname(f)
        file_name = os.path.basename(f)
        bknd, bknd_name, algebra, setup = __read_context_metadata(
            logging, input_dir, file_name
        )

        if bknd not in benchmarks:
            benchmarks[bknd] = {}

        if bknd_name not in benchmarks[bknd]:
            benchmarks[bknd][bknd_name] = []

        context, _ = read_benchmark_data(logging, input_dir, file_name)
        if bench_type in context["executable"]:
            benchmarks[bknd][bknd_name].append(
                benchmark_metadata(name=bknd_name, algebra=algebra, setup=setup, file=f)
            )

    # Register benchmarks to be run
    for bknd, metadata_dict in benchmarks.items():
        # No benchmarks to be run
        if len(algebra_plugins) == 0:
            break

        # Try to find the processor name
        bknd_name = "Unknown"
        if bknd == "CUDA" or bknd == "SYCL":
            # Try to get the GPU name
            try:
                gpu_str = str(
                    subprocess.check_output(
                        ["nvidia-smi", "--query-gpu", "name", "--format=csv,noheader"]
                    )
                )
                gpu_str = __compactify_bknd_name(gpu_str[2:])

                # Strip some unwanted characters
                if gpu_str[-1] == '"' or gpu_str[-1] == "'":
                    gpu_str = gpu_str[:-1]
                if gpu_str[-2:] == "\\n":
                    gpu_str = gpu_str[:-2]

                if len(gpu_str) != 0:
                    bknd_name = gpu_str.rstrip(os.linesep)

            except Exception as e:
                # Name remains 'Unknown'
                print(e)
        else:
            bknd_name = __compactify_bknd_name(platform.processor())

        if bknd_name not in metadata_dict:
            metadata_dict[bknd_name] = []

        # The algebra plugins that have already been registered from input file
        registered_algebra = [data.algebra for data in metadata_dict[bknd_name]]

        # The requested algebra-plugins
        for algebra in algebra_plugins:
            # Prefer reading from file if it has been provided by the user
            if algebra in registered_algebra:
                continue

            # Parse the detector setup
            setup = ""
            add_delim = lambda s: s + ", "
            if not args.grid_file:
                setup = setup + "no grids"
            if not args.material_file:
                if len(setup) != 0:
                    setup = add_delim(setup)
                setup = setup + "no mat."
            if not args.covariance_transport:
                if len(setup) != 0:
                    setup = add_delim(setup)
                setup = setup + "no cov."

            binary = (
                f"{bindir}/detray_propagation_{bench_type}_{bknd.lower()}_{algebra}"
            )
            file_bknd_name = bknd_name.replace(" ", "_")
            file_setup = setup.replace(" ", "_")
            file_setup = file_setup.replace(",", "")
            file_setup = file_setup.replace(".", "")
            data_file = f"{det_name}_{bench_type}_{bknd.lower()}_{file_bknd_name}_{algebra}_{file_setup}.json"

            metadata = benchmark_metadata(
                name=bknd_name, algebra=algebra, setup=setup, file=data_file
            )

            # If the results should not be read from file, run the benchmark
            if data_file not in (os.path.basename(f) for f in input_data_files):
                # Register binary if it exists
                if os.path.isdir(bindir) and os.path.isfile(binary):
                    metadata = metadata._replace(bin=binary)
                else:
                    logging.warning(
                        f"Propagation {bench_type} binary not found! ({binary})"
                    )
                    continue

            metadata_dict[bknd_name].append(metadata)

    return benchmarks


# Run all benchmarks in 'benchmark_dict' that were registered with a binary file
def __run_benchmarks(benchmark_dict, args_list, benchmark_options):
    for bknd, metadata_dict in benchmark_dict.items():
        for bknd_name, metadata_list in metadata_dict.items():
            for metadata in metadata_list:
                if metadata.bin is not None:
                    subprocess.run(
                        [
                            metadata.bin,
                            f"--bknd_name={bknd_name}",
                            f"--benchmark_out=./{metadata.file}",
                        ]
                        + benchmark_options
                        + args_list
                    )


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
        default=[],
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

    # Check user provided benchmark result files
    input_data_files = __parse_input_data_files(args, logging)

    # Unique set of algebra plugins to be included in the plots
    algebra_plugins = set(args.algebra_plugins)

    if len(algebra_plugins) == 0 and len(input_data_files) == 0:
        logging.error(
            "No data file for plotting and no algebra plugins specified to run benchmarks for. Quitting"
        )
        sys.exit(1)

    # Get dictionary of benchmark files per hardware backend type
    benchmarks = __generate_benchmark_dict(
        args, logging, bindir, det_name, input_data_files, algebra_plugins, "benchmark"
    )

    scaling = __generate_benchmark_dict(
        args, logging, bindir, det_name, input_data_files, algebra_plugins, "scaling"
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
    ]

    # Run requested benchmark and scaling tests
    __run_benchmarks(benchmarks, args_list, benchmark_options)
    __run_benchmarks(scaling, args_list, benchmark_options)

    # ----------------------------------------------------------------------plot

    logging.info("Generating plots...\n")

    plot_factory = plt_factory(out_dir, logging)

    def make_label(algebra, setup):
        label = f"{algebra}"
        if len(setup) != 0:
            label = label + f" ({setup})"
        return label

    # Plot all data files per hardware backend
    # (comparison of different algebra-plugins)
    for bknd, metadata_dict in benchmarks.items():
        for bknd_name, metadata_list in metadata_dict.items():
            for metadata in metadata_list:
                # Get file list and plot labels
                files = [metadata.file for metadata in metadata_list]
                plot_labels = [
                    make_label(metadata.algebra, metadata.setup)
                    for metadata in metadata_list
                ]

                file_bknd_name = bknd_name.replace(" ", "_")
                plot_benchmark_data(
                    logging,
                    input_dir,
                    det_name,
                    files,
                    plot_labels,
                    f"hardware backend: {bknd} ({bknd_name})",
                    f"prop_benchmark_algebra-plugin_comparison_{bknd}_{file_bknd_name}",
                    plot_factory,
                    out_format,
                )

    # Plot results for different hardware backends using the same algebra plugin
    # (comparison of different hardware backends)
    discovered_algs = []
    for bknd, metadata_dict in benchmarks.items():
        for bknd_name, metadata_list in metadata_dict.items():
            for metadata in metadata_list:
                discovered_algs.append(metadata.algebra)

    discovered_algs = set(discovered_algs)

    # Build up the list of data files and corresponding plot labels
    for algebra in discovered_algs:
        data_files_per_plugin = []
        plot_labels = []

        for bknd, metadata_dict in benchmarks.items():

            for bknd_name, metadata_list in metadata_dict.items():

                for metadata in metadata_list:
                    if algebra == metadata.algebra:
                        data_files_per_plugin.append(metadata.file)
                        label = make_label(f"{bknd}: {bknd_name}", metadata.setup)
                        plot_labels.append(label)

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

    # Plot all data files per hardware backend for propagation scaling
    for bknd, metadata_dict in scaling.items():
        for bknd_name, metadata_list in metadata_dict.items():
            for metadata in metadata_list:
                # Get file list and plot labels
                files = [metadata.file for metadata in metadata_list]
                plot_labels = [
                    make_label(metadata.algebra, metadata.setup)
                    for metadata in metadata_list
                ]

                plot_scaling_data(
                    logging,
                    input_dir,
                    det_name,
                    files,
                    plot_labels,
                    f"hardware backend: {bknd} ({bknd_name})",
                    plot_factory,
                    out_format,
                )


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
