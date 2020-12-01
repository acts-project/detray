add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

set (BENCHMARK_ENABLE_INSTALL OFF)
set (BENCHMARK_ENABLE_TESTING OFF)
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/benchmark" "extern/benchmark")

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

include(GoogleTest)

macro(add_detray_test TESTNAME FILES LIBRARIES)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main detray::core ${LIBRARIES})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    add_test(unit_test_${TESTNAME} ${TESTNAME})
    set_tests_properties(unit_test_${TESTNAME} PROPERTIES 
        ENVIRONMENT DETRAY_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/tests/data)
endmacro()

macro(add_detray_benchmark TESTNAME FILES LIBRARIES)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} benchmark detray::core ${LIBRARIES})
    #add_test(benchmark_${TESTNAME} ${TESTNAME})
endmacro()

add_subdirectory(core)
if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()