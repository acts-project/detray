# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the SYCL build flags.
include(detray-compiler-options-sycl)

# Enable SYCL as a language.
enable_language(SYCL)

# Set up the library.
#detray_add_library( detray_sycl_test sycl_test TYPE SHARED
#    "include/detray/utils/calculate1DimNdRange.hpp"
#    "include/detray/utils/get_queue.hpp"
#    "include/detray/utils/queue_wrapper.hpp"
#    "src/utils/calculate1DimNdRange.sycl"
    #"src/utils/get_queue.hpp"
    #"src/utils/get_queue.sycl"
#)

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebras "eigen")
 endif()

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
 