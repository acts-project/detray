# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Set the CPU build flags.
include(detray-compiler-options-cpp)

# Set up a test library, which the "new style" benchmarks and tests could use.
add_library(detray_benchmark_cpu INTERFACE "propagation_benchmark.hpp")

add_library(detray::benchmark_cpu ALIAS detray_benchmark_cpu)

target_link_libraries(detray_benchmark_cpu INTERFACE detray::benchmarks)
