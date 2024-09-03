# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Set up a benchamrk library for CUDA
add_library( detray_benchmark_cuda STATIC "propagation_benchmark.hpp"
                                          "propagation_benchmark.cu" )

add_library( detray::benchmark_cuda ALIAS detray_benchmark_cuda )

target_link_libraries( detray_benchmark_cuda PUBLIC detray::benchmarks detray::core_array detray::test_common vecmem::cuda detray::utils )
