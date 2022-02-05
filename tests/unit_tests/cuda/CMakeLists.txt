# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable CUDA as a language.
enable_language( CUDA )

# Set the CUDA build flags.
include( detray-compiler-options-cuda )

# "Core" tests.
detray_add_test( core_cuda
   "grids_grid2_cuda.cpp" "grids_grid2_cuda_kernel.hpp"
   "grids_grid2_cuda_kernel.cu" "utils_enumerate_cuda.cpp"
   "utils_enumerate_cuda_kernel.hpp" "utils_enumerate_cuda_kernel.cu"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common vecmem::cuda )

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set( algebras "array" )
if( DETRAY_EIGEN_PLUGIN )
   list( APPEND algebras "eigen" )
endif()

foreach(algebra ${algebras})

   # Unit tests for the selected algebra.
   detray_add_test( ${algebra}_cuda
      "transform_store_cuda.cpp" "transform_store_cuda_kernel.hpp"
      "transform_store_cuda_kernel.cu"
      "mask_store_cuda.cpp" "mask_store_cuda_kernel.hpp"
      "mask_store_cuda_kernel.cu"
      "tuple_test_cuda.cpp" "tuple_test_cuda_kernel.hpp"
      "tuple_test_cuda_kernel.cu"
      "surfaces_finder_cuda.cpp" "surfaces_finder_cuda_kernel.hpp"
      "surfaces_finder_cuda_kernel.cu"
      "detector_cuda.cpp" "detector_cuda_kernel.hpp"
      "detector_cuda_kernel.cu"
      "navigator_cuda.cpp" "navigator_cuda_kernel.hpp"
      "navigator_cuda_kernel.cu"
      LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray_tests_common
                     detray::${algebra} )
   target_compile_definitions( detray_test_${algebra}_cuda
      PRIVATE ${algebra}=${algebra} )

endforeach()
