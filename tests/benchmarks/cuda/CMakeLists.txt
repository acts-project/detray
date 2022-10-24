# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required( VERSION 3.18 )

# Enable CUDA as a language.
enable_language( CUDA )

# Set the CUDA build flags.
include( detray-compiler-options-cuda )

# Look for openMP, which is used for the CPU benchmark
find_package(OpenMP)

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set( algebras "array" )
if( DETRAY_EIGEN_PLUGIN )
   list( APPEND algebras "eigen" )
endif()

foreach(algebra ${algebras})

detray_add_executable( ${algebra}_propagator_cuda
   "benchmark_propagator_cuda.cpp" "benchmark_propagator_cuda_kernel.hpp"
   "benchmark_propagator_cuda_kernel.cu" 
   LINK_LIBRARIES benchmark::benchmark detray_tests_common detray::${algebra} vecmem::cuda OpenMP::OpenMP_CXX)

   target_compile_definitions( detray_${algebra}_propagator_cuda
      PRIVATE ${algebra}=${algebra} )

endforeach()
