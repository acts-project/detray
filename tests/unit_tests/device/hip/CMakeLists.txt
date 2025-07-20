# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message (STATUS "Building detray HIP unit tests")

cmake_minimum_required(VERSION 3.21) # HIP langauge support requires minimum 3.21 
#cmake_policy(VERSION 3.21.3...3.27)
#find_package(hip REQUIRED)

# Enable HIP as a language.
enable_language(HIP)


# Set the HIP build flags.
include(detray-compiler-options-hip)
set_source_files_properties(detector_hip.cpp PROPERTIES LANGUAGE HIP)

set_source_files_properties(detector_hip_kernel.cpp PROPERTIES LANGUAGE HIP)

add_executable(detray_build
    detector_hip.cpp
    detector_hip_kernel.cpp
)

target_sources(detray_build PRIVATE
    detector_hip.cpp
    detector_hip_kernel.cpp
)


#target_include_directories(detray_build PRIVATE
#  ${CMAKE_CURRENT_SOURCE_DIR}
#  ${HIP_PATH}/include
#)

# Define the HIP platform
target_compile_definitions(detray_build PRIVATE __HIP_PLATFORM_AMD__)


# Define the algebra plugin explicitly (ex: array)
target_compile_definitions(detray_build PRIVATE DETRAY_USE_ALGEBRA_ARRAY)

# Link with the algebra plugin and other required libraries
target_link_libraries(detray_build
    PRIVATE
        vecmem::hip
        GTest::gtest_main
        detray::core
        detray::algebra_array
        detray::test_common
        detray::test_utils
)

# Add project headers


# Optional: include HIP headers via module
if(DEFINED ENV{HIP_PATH})
    target_include_directories(detray_build PRIVATE $ENV{HIP_PATH}/include)
endif()



# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
set(algebras "array")

if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebras "eigen")
endif()

foreach(algebra ${algebras})
    # Unit tests for the selected algebra.
    #detray_add_unit_test(hip_${algebra}
    #   "detector_hip_kernel.hpp"
    #   "detector_hip.cpp"
    #   "detector_hip_kernel.cpp"
    detray_add_unit_test(hip_${algebra}
       "detector_hip_kernel.hpp"
       "detector_hip.cpp"
       "detector_hip_kernel.cpp"
       LINK_LIBRARIES GTest::gtest_main vecmem::hip  detray::core 
       detray::algebra_${algebra} detray::test_common detray::test_utils
    )

    # Ajout des définitions pour HIP + algebra
    target_compile_definitions(
        detray_unit_test_hip_${algebra}
        PRIVATE ${algebra}=${algebra} __HIP_PLATFORM_AMD__
    )



endforeach()