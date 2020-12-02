set (BENCHMARK_ENABLE_INSTALL OFF)
set (BENCHMARK_ENABLE_TESTING OFF)
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/benchmark" "extern/benchmark")

macro(add_detray_benchmark TESTNAME FILES PLUGIN_LIBRARY)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} benchmark detray::core detray::tests_common detray::tests_io ${PLUGIN_LIBRARY})
    add_test(benchmark_${TESTNAME} ${TESTNAME})
    set_tests_properties(benchmark_${TESTNAME} PROPERTIES 
        ENVIRONMENT DETRAY_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/tests/data)
endmacro()

if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()