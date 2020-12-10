add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

include(GoogleTest)

macro(add_detray_test TESTNAME FILES PLUGIN_LIBRARY)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    target_link_libraries(${TESTNAME} detray::core detray::tests_common detray::tests_io)
    target_link_libraries(${TESTNAME} ${PLUGIN_LIBRARY})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    add_test(unit_test_${TESTNAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TESTNAME})
    set_tests_properties(unit_test_${TESTNAME} PROPERTIES 
        ENVIRONMENT DETRAY_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/tests/data)
endmacro()

add_subdirectory(core)
if(DETRAY_ARRAY_PLUGIN)
    add_subdirectory(array)
endif()
if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()
