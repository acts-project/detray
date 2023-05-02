# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up a helper library, which would be used by most of the tests.
add_library( detray_tests_common INTERFACE )
target_include_directories( detray_tests_common
   INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include" )
target_link_libraries( detray_tests_common
   INTERFACE vecmem::core detray::core detray::io detray::utils covfie::core )

if( DETRAY_BENCHMARKS AND DETRAY_BENCHMARKS_REP )
   message( STATUS "Repetitions for benchmarks: ${DETRAY_BENCHMARKS_REP}" )
   target_compile_definitions( detray_tests_common
      INTERFACE DETRAY_BENCHMARKS_REP=${DETRAY_BENCHMARKS_REP} )
endif()

if( DETRAY_BENCHMARKS_MULTITHREAD )
   message( STATUS "Using multithreaded benchmarks" )
   target_compile_definitions( detray_tests_common
      INTERFACE DETRAY_BENCHMARKS_MULTITHREAD )
endif()
