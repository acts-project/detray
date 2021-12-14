# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::array' plugin")

# Set up the library.
detray_add_library( detray_array array
   "include/detray/plugins/algebra/array_definitions.hpp" )
target_link_libraries( detray_array
   INTERFACE algebra::array_cmath vecmem::core )
target_compile_definitions( detray_array
   INTERFACE DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE} )
