# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

#----------------------------- no errrors
message(STATUS "Building detray HIP integration tests")

#HIP langauge support requires minimum 3.21
cmake_minimum_required(VERSION 3.21)

#set(CMAKE_HIP_COMPILER_FORCED TRUE)
#set(CMAKE_HIP_STANDARD_COMPUTED_DEFAULT 17)
# Enable HIP as a language.
enable_language(HIP)

# Set the HIP build flags.
include(detray-compiler-options-hip)

#---------------------------------------------

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebras "eigen")
endif()

# error in covfie::hip : No known features for HIP compiler "AMD" Version 6.....

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    detray_add_integration_test(hip_${algebra}
       "propagator_hip_kernel.hpp"
       "propagator_hip.hip"
       "propagator_hip_kernel.hip"
       LINK_LIBRARIES GTest::gtest_main vecmem::hip #covfie::hip
                       detray::core detray::algebra_${algebra} detray::test_hip
                       detray::test_common
    )

    target_compile_definitions(
        detray_integration_test_hip_${algebra}
        PRIVATE ${algebra}=${algebra}
    )

    set_tests_properties(
        detray_integration_test_hip_${algebra}
        PROPERTIES
            DEPENDS "detray_unit_test_hip;detray_unit_test_hip_${algebra}"
    )
endforeach()
