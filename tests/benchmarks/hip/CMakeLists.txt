# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray HIP benchmarks")

cmake_minimum_required(VERSION 3.21) # HIP langauge support requires minimum 3.21

find_package(HIPToolkit)

# Enable HIP as a language.
enable_language(HIP)

# Set the HIP build flags.
include(detray-compiler-options-hip)

# Build benchmarks for multiple algebra plugins
# Currently vc and smatrix is not supported on device
set(algebra_plugins "array")
if(DETRAY_EIGEN_PLUGIN)
    list(APPEND algebra_plugins "eigen")
endif()

foreach(algebra ${algebra_plugins})
    detray_add_executable(benchmark_hip_propagation_${CMAKE_HIP_PLATFORM}_${algebra}
      "propagation.cpp"
       LINK_LIBRARIES detray::benchmark_hip_${CMAKE_HIP_PLATFORM}_${algebra} detray::core_${algebra} vecmem::hip detray::test_common HIP::hiprt
    )

    if(
        ("${CMAKE_HIP_PLATFORM}" STREQUAL "hcc")
        OR ("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
    )
        set_source_files_properties(propagation.cpp PROPERTIES LANGUAGE HIP)
    endif()

    target_compile_definitions(
        detray_benchmark_hip_propagation_${CMAKE_HIP_PLATFORM}_${algebra}
        PRIVATE ${algebra}=${algebra}
    )

    target_compile_options(
        detray_benchmark_hip_propagation_${CMAKE_HIP_PLATFORM}_${algebra}
        PRIVATE "-march=native" "-ftree-vectorize"
    )
endforeach()
