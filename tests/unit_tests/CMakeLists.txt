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
  if ( ARG_LINK_LIBRARIES )
    target_link_libraries( ${test_exe_name} PRIVATE ${ARG_LINK_LIBRARIES})
    target_link_libraries( ${test_exe_name} PRIVATE detray::core detray::io detray::tests_common)
    target_link_libraries( ${test_exe_name} PRIVATE gtest gmock gtest_main)
  endif()
  foreach( _config "" "_DEBUG" "_RELEASE" "_MINSIZEREL" "_RELWITHDEBINFO" )
    set_property( TARGET ${test_exe_name} PROPERTY
      RUNTIME_OUTPUT_DIRECTORY${_config} "${CMAKE_BINARY_DIR}/test-bin" )
  endforeach()

  # Check if it is cuda test
  if(${test_exe_name} MATCHES ".*_cuda")
    target_compile_options(${test_exe_name}
      PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
      --expt-relaxed-constexpr>)
  endif()
  
  # Run the executable as the test.
  add_test( NAME ${test_exe_name}
    COMMAND $<TARGET_FILE:${test_exe_name}> )

endfunction(add_detray_test)  
  
set(all_unit_tests "core")
list(APPEND all_unit_tests "annulus2;cylinder3;rectangle2;ring2;single3;unmasked;trapezoid2")
list(APPEND all_unit_tests "propagator;line_stepper;navigator;cylinder_intersection;planar_intersection;intersection_kernel")
list(APPEND all_unit_tests "detector;transform_store;masks_container")

add_subdirectory(core)
if(DETRAY_ARRAY_PLUGIN)
    add_subdirectory(array)
endif()
if(DETRAY_EIGEN_PLUGIN)
    add_subdirectory(eigen)
endif()
if(DETRAY_SMATRIX_PLUGIN)
    add_subdirectory(smatrix)
endif()
if(DETRAY_VC_ALGEBRA_PLUGIN)
    add_subdirectory(vc)
endif()
if(DETRAY_BUILD_CUDA)
  add_subdirectory(cuda)
endif()
