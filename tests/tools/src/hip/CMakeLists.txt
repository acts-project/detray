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



if(DETRAY_BUILD_BENCHMARKS)
    # Build benchmarks for multiple algebra plugins
    # Currently vc and smatrix is not supported on device
    set(algebra_plugins "array")
    if(DETRAY_EIGEN_PLUGIN)
        list(APPEND algebra_plugins "eigen")
    endif()

    foreach(algebra ${algebra_plugins})
        detray_add_executable(propagation_benchmark_hip-${CMAKE_HIP_PLATFORM}_${algebra}
        "propagation_benchmark_hip.cpp"
        LINK_LIBRARIES detray::benchmark_hip_${algebra} vecmem::hip detray::tools detray::test_common HIP::hiprt
        )
    endforeach()
endif()
