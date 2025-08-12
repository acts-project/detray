# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray host benchmarks")

# Set up the benchmarking options.
option(DETRAY_BENCHMARK_MULTITHREAD "Enable multithreaded benchmarks" OFF)
option(DETRAY_BENCHMARK_PRINTOUTS "Enable printouts in the benchmarks" OFF)

# Look for openMP, which is used for the CPU propagation benchmark
find_package(OpenMP)

# Macro setting up the CPU benchmarks for a specific algebra plugin.
macro(detray_add_cpu_benchmark algebra)
    # Build the benchmark executable.
    detray_add_executable(benchmark_cpu_${algebra}
       "find_volume.cpp"
       "grid.cpp"
       "grid2.cpp"
       "intersect_all.cpp"
       "intersect_surfaces.cpp"
       "masks.cpp"
       LINK_LIBRARIES benchmark::benchmark benchmark::benchmark_main vecmem::core detray::benchmarks
                      detray::core_${algebra} detray::test_common
    )

    target_compile_options(
        detray_benchmark_cpu_${algebra}
        PRIVATE "-march=native"
    )

    # Set the benchmark specific compilation options.
    if(DETRAY_BENCHMARKS_MULTITHREAD)
        target_compile_definitions(
            detray_benchmark_cpu_${algebra}
            PRIVATE DETRAY_BENCHMARKS_MULTITHREAD
        )
    endif()

    if(DETRAY_BENCHMARK_PRINTOUTS)
        target_compile_definitions(
            detray_benchmark_cpu_${algebra}
            PRIVATE DETRAY_BENCHMARK_PRINTOUTS
        )
    endif()

    # Build the benchmark executable for the propagation
    detray_add_executable( benchmark_cpu_propagation_${algebra}
       "propagation.cpp"
       LINK_LIBRARIES detray::benchmark_cpu benchmark::benchmark_main
                     vecmem::core detray::core_${algebra} detray::test_common
    )

    target_compile_options(
        detray_benchmark_cpu_propagation_${algebra}
        PRIVATE "-march=native" "-ftree-vectorize"
    )

    if(OpenMP_CXX_FOUND)
        target_link_libraries(
            detray_benchmark_cpu_propagation_${algebra}
            PRIVATE OpenMP::OpenMP_CXX
        )
    endif()
endmacro()

# Build the array benchmark.
detray_add_cpu_benchmark( array )

# Build the Eigen benchmark executable.
if(DETRAY_EIGEN_PLUGIN)
    detray_add_cpu_benchmark( eigen )
endif()

# Build the SMatrix benchmark executable.
if(DETRAY_SMATRIX_PLUGIN)
    detray_add_cpu_benchmark( smatrix )
endif()

# Build the Vc benchmark executable.
if(DETRAY_VC_AOS_PLUGIN)
    detray_add_cpu_benchmark( vc_aos )
endif()

if(DETRAY_VC_SOA_PLUGIN)
    macro(detray_add_soa_benchmark algebra)
        # Build the benchmark executable.
        detray_add_executable(benchmark_cpu_vc_soa_vs_${algebra}
         "intersectors.cpp"
           LINK_LIBRARIES detray::benchmark_cpu benchmark::benchmark benchmark::benchmark_main vecmem::core detray::core_vc_soa detray::core_vc_aos detray::core_${algebra}
           detray::test_common
        )

        target_compile_options(
            detray_benchmark_cpu_vc_soa_vs_${algebra}
            PRIVATE "-march=native"
        )

        # Add Vc specific flags
        target_compile_options(
            detray_benchmark_cpu_vc_soa_vs_${algebra}
            INTERFACE ${Vc_COMPILE_FLAGS} ${Vc_ARCHITECTURE_FLAGS}
        )

        # Set the benchmark specific compilation options.
        if(DETRAY_BENCHMARKS_MULTITHREAD)
            target_compile_definitions(
                detray_benchmark_cpu_vc_soa_vs_${algebra}
                PRIVATE DETRAY_BENCHMARKS_MULTITHREAD
            )
        endif()

        if(DETRAY_BENCHMARK_PRINTOUTS)
            target_compile_definitions(
                detray_benchmark_cpu_vc_soa_vs_${algebra}
                PRIVATE DETRAY_BENCHMARK_PRINTOUTS
            )
        endif()
    endmacro()

    detray_add_soa_benchmark( array )
    # detray_add_soa_benchmark( eigen )
    # detray_add_soa_benchmark( vc_aos )
endif()
