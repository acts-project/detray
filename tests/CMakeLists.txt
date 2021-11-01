
add_subdirectory(common)
add_subdirectory(unit_tests)

if(DETRAY_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
