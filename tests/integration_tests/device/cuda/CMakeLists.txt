# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# C++17 support for CUDA requires CMake 3.18.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)

# Set the CUDA build flags.
include(detray-compiler-options-cuda)

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
   list(APPEND algebras "eigen")
endif()

foreach(algebra ${algebras})
   # Unit tests for the selected algebra.
   detray_add_integration_test(cuda_${algebra}
      "propagator_cuda_kernel.hpp"
      "propagator_cuda.cpp"
      "propagator_cuda_kernel.cu"
      LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray::test covfie::cuda
                     detray::core detray::algebra_${algebra} detray::utils )

   target_compile_definitions(detray_integration_test_cuda_${algebra}
      PRIVATE ${algebra}=${algebra})

   set_tests_properties(detray_integration_test_cuda_${algebra} 
                        PROPERTIES DEPENDS
                     "detray_unit_test_cuda;detray_unit_test_cuda_${algebra}")
endforeach()
