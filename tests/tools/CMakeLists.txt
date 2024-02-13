# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

include( CMakeFindDependencyMacro )

find_dependency( Boost COMPONENTS program_options )

add_library( detray_tools INTERFACE )
add_library( detray::tools ALIAS detray_tools )

target_include_directories( detray_tools
   INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include" )

target_link_libraries( detray_tools
   INTERFACE vecmem::core detray::test detray::io detray::core_array )

add_subdirectory( src )
