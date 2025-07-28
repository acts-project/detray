# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0


message (STATUS "Building detray HIP unit tests")
cmake_minimum_required(VERSION 3.21) # HIP langauge support requires minimum 3.21 

# Enable HIP as a language.
enable_language(HIP)


# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebras "eigen")
endif()
foreach(algebra ${algebras})
    detray_add_unit_test(hip_${algebra}
       "detector_hip_kernel.hpp"
       "detector_hip.hip"
       "detector_hip_kernel.hip"
       
       LINK_LIBRARIES GTest::gtest_main vecmem::hip  detray::core 
       detray::algebra_${algebra} detray::test_common detray::test_utils
    )

    # Add definitions for HIP + algebra
    target_compile_definitions(
        detray_unit_test_hip_${algebra}
        PRIVATE ${algebra}=${algebra}
    )



endforeach()