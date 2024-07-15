# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the CPU build flags.
include(detray-compiler-options-cpp)

# Set up a benchamrk library for CPU
add_library( detray_benchmark_cpu INTERFACE "propagation_benchmark.hpp")

add_library( detray::benchmark_cpu ALIAS detray_benchmark_cpu )

target_link_libraries( detray_benchmark_cpu INTERFACE detray::benchmarks detray::test_common )
