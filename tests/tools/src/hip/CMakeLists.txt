#Detray library, part of the ACTS project(R& D line)
#
#(c) 2024-2025 CERN for the benefit of the ACTS project
#
#Mozilla Public License Version 2.0

message(STATUS "Building detray device command line tools")

find_package(HIPToolkit)

# Enable HIP as a language.
enable_language(HIP)

include(CMakeFindDependencyMacro)


find_dependency(Boost COMPONENTS program_options REQUIRED)

if(DETRAY_BUILD_TESTING)
    # Build the navigation validation executable.
    detray_add_executable(navigation_validation_hip_${CMAKE_HIP_PLATFORM}
                        "navigation_validation_hip.cpp"
                        LINK_LIBRARIES HIP::hiprt  GTest::gtest GTest::gtest_main
                        Boost::program_options detray::test_hip_${CMAKE_HIP_PLATFORM} detray::tools
    )
    if(
        ("${CMAKE_HIP_PLATFORM}" STREQUAL "hcc")
        OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
    )
        set_source_files_properties(navigation_validation_hip.cpp PROPERTIES LANGUAGE HIP)
    endif()

    # Build the material hip validation executable.
    detray_add_executable(material_validation_hip_${CMAKE_HIP_PLATFORM}
                        "material_validation_hip.cpp"
                        LINK_LIBRARIES HIP::hiprt  GTest::gtest GTest::gtest_main
                        Boost::program_options detray::test_hip_${CMAKE_HIP_PLATFORM} detray::tools
    )
    if(
        ("${CMAKE_HIP_PLATFORM}" STREQUAL "hcc")
        OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
    )
        set_source_files_properties(material_validation_hip.cpp PROPERTIES LANGUAGE HIP)
    endif()
endif()

if(DETRAY_BUILD_BENCHMARKS)
    # Build benchmarks for multiple algebra plugins
    # Currently vc and smatrix is not supported on device
    set(algebra_plugins "array")
    if(DETRAY_EIGEN_PLUGIN)
        list(APPEND algebra_plugins "eigen")
    endif()

    foreach(algebra ${algebra_plugins})
        detray_add_executable(propagation_benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra}
        "propagation_benchmark_hip.cpp"
        LINK_LIBRARIES detray::benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra} HIP::hiprt vecmem::hip detray::tools detray::test_common 
        )
        # Make hipcc interpret .cpp files as .hip files to make linking work
        # (only needed for the amd backend)
        if(
            ("${CMAKE_HIP_PLATFORM}" STREQUAL "hcc")
            OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
        )
            set_source_files_properties(propagation_benchmark_hip.cpp PROPERTIES LANGUAGE HIP)
        endif()
        target_compile_definitions(
            detray_propagation_benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra}

            PRIVATE HIP_PLATFORM_NAME="${CMAKE_HIP_PLATFORM}"
        )

    endforeach()
endif()
