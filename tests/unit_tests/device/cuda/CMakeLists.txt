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

# "Core" tests.
detray_add_test(core_cuda
   "single_store.cu"
   "container_cuda.cpp" "container_cuda_kernel.hpp" "container_cuda_kernel.cu"
   "grids_grid2_cuda.cpp" "grids_grid2_cuda_kernel.hpp"
   "grids_grid2_cuda_kernel.cu" "utils_ranges_cuda.cpp"
   "utils_ranges_cuda_kernel.hpp" "utils_ranges_cuda_kernel.cu"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common vecmem::cuda
                  detray::core_array )

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
   list(APPEND algebras "eigen")
endif()

foreach(algebra ${algebras})
   # Unit tests for the selected algebra.
   detray_add_test(${algebra}_cuda
      "transform_store_cuda.cpp" "transform_store_cuda_kernel.hpp"
      "transform_store_cuda_kernel.cu"
      "mask_store_cuda.cpp" "mask_store_cuda_kernel.hpp"
      "mask_store_cuda_kernel.cu"
      "detector_cuda.cpp" "detector_cuda_kernel.hpp"
      "detector_cuda_kernel.cu"
      "navigator_cuda.cpp" "navigator_cuda_kernel.hpp"
      "navigator_cuda_kernel.cu"
      "sf_finders_grid_cuda.cpp" "sf_finders_grid_cuda_kernel.hpp"
      "sf_finders_grid_cuda_kernel.cu"
      "propagator_cuda.cpp" "propagator_cuda_kernel.hpp"
      "propagator_cuda_kernel.cu"
      LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray_tests_common
      detray::algebra_${algebra})
   target_compile_definitions(detray_test_${algebra}_cuda
      PRIVATE ${algebra}=${algebra})
endforeach()
