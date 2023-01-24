# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the SYCL build flags.
include(detray-compiler-options-sycl)

# Enable SYCL as a language.
enable_language(SYCL)

# make unit tests for multiple algebras
# Currently eigen, vc and smatrix is not supported
set(algebras "array")

foreach(algebra ${algebras})
   # Unit tests for the selected algebra.
   detray_add_test(${algebra}_sycl
      "propagator.sycl" "propagator_kernel.hpp"
      "propagator_kernel.sycl"
      LINK_LIBRARIES GTest::gtest_main vecmem::sycl detray_tests_common
      detray::${algebra})
   target_compile_definitions(detray_test_${algebra}_sycl
      PRIVATE ${algebra}=${algebra})
endforeach()
