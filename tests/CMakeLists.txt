# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the common C++ flags.
include( detray-compiler-options-cpp )
include_directories( SYSTEM $<TARGET_PROPERTY:covfie::core,INTERFACE_INCLUDE_DIRECTORIES> )

# Set up a test library, which the "new style" benchmarks and tests could use.
file( GLOB _detray_test_headers
   RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "include/detray/test/*.hpp"
   "include/detray/test/detail/*.hpp" 
   "include/detray/test/intersection/*.hpp"
   "include/detray/test/utils/*.hpp")

add_library( detray_test INTERFACE "${_detray_test_headers}" )
add_library( detray::test ALIAS detray_test )

target_include_directories( detray_test
   INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include" )

target_link_libraries( detray_test INTERFACE vecmem::core detray::core )

unset( _detray_test_headers )

# Include all of the code-holding sub-directories.
if( DETRAY_BENCHMARKS )
   add_subdirectory( benchmarks )
endif()
add_subdirectory( integration_tests )
add_subdirectory( tools )
add_subdirectory( unit_tests )
