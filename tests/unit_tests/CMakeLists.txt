message(STATUS "Testing: 'detray::unit_tests' enabled")

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

include(GoogleTest)

function(add_detray_test name)
  
  # Parse the function's options.
  cmake_parse_arguments( ARG "" "" "LINK_LIBRARIES" ${ARGN} )

  # Create the test executable.
  set( test_exe_name "detray_test_${name}" )
  add_executable(${test_exe_name} ${ARG_UNPARSED_ARGUMENTS})
  target_link_libraries( ${test_exe_name} PRIVATE gtest gmock gtest_main)
  target_link_libraries( ${test_exe_name} PRIVATE detray::core detray::io detray::tests_common)  
  if ( ARG_LINK_LIBRARIES )
    target_link_libraries( ${test_exe_name} PRIVATE ${ARG_LINK_LIBRARIES})
  endif()
  foreach( _config "" "_DEBUG" "_RELEASE" "_MINSIZEREL" "_RELWITHDEBINFO" )
    set_property( TARGET ${test_exe_name} PROPERTY
      RUNTIME_OUTPUT_DIRECTORY${_config} "${CMAKE_BINARY_DIR}/test-bin" )
  endforeach()

  target_compile_options(${test_exe_name}
    PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
    --expt-relaxed-constexpr>)
  
  # Run the executable as the test.
  add_test( NAME ${test_exe_name} COMMAND $<TARGET_FILE:${test_exe_name}> )

  set_tests_properties(${test_exe_name} PROPERTIES ENVIRONMENT DETRAY_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data/)
  
endfunction(add_detray_test)  
  
set(all_unit_tests "core")
list(APPEND all_unit_tests "annulus2;cylinder3;rectangle2;ring2;single3;unmasked;trapezoid2")
list(APPEND all_unit_tests "propagator;line_stepper;navigator;single_type_navigator;cylinder_intersection;planar_intersection;intersection_kernel")
list(APPEND all_unit_tests "detector;transform_store;masks_container;mask_store")
list(APPEND all_unit_tests "volume;index_geometry;unified_index_geometry;geometry_scan;geometry_navigation;toy_geometry;geometry_linking")

add_subdirectory(core)
add_subdirectory(array)
if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()
if(DETRAY_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()
if(DETRAY_VC_PLUGIN)
    add_subdirectory(vc)
endif()
if(DETRAY_BUILD_CUDA)
  add_subdirectory(cuda)
endif()
