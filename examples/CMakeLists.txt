# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the common C++ flags.
include( detray-compiler-options-cpp )
include_directories( SYSTEM $<TARGET_PROPERTY:dfelibs::dfelibs,INTERFACE_INCLUDE_DIRECTORIES> )

# Set up a example library, which showcases how to use detray components.
file( GLOB _detray_example_headers
   RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "include/detray/examples/*.hpp" )
add_library( detray_example INTERFACE "${_detray_example_headers}" )
add_library( detray::example ALIAS detray_example )
target_include_directories( detray_example
   INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include" )
target_link_libraries( detray_example INTERFACE
   vecmem::core detray::core )
unset( _detray_example_headers )

# Include all of the code-holding sub-directories.
add_subdirectory( src )
