# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# Examples using the array algebra plugin.
detray_add_tutorial(propagation_cuda
   "propagation.cpp"
   "propagation.hpp"
   "propagation_kernel.cu"
   LINK_LIBRARIES vecmem::cuda detray::core_array covfie::core
                  detray::tutorial detray::utils)
target_compile_definitions(detray_tutorial_propagation_cuda PRIVATE array=array)
