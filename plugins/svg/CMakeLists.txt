# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::svg' component")

include_directories( SYSTEM $<TARGET_PROPERTY:actsvg::meta,INTERFACE_INCLUDE_DIRECTORIES> )
include_directories( SYSTEM $<TARGET_PROPERTY:actsvg::core,INTERFACE_INCLUDE_DIRECTORIES> )

# Set up the library.
file( GLOB _detray_svg_public_headers
   RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "include/detray/plugins/actsvg_visualization/*.hpp" )

detray_add_library( detray_svg svg
   ${_detray_svg_public_headers} )
target_link_libraries( detray_svg
   INTERFACE detray::core actsvg::core actsvg::meta )

# Clean up.
unset( _detray_svg_public_headers )