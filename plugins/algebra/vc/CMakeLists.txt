# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::vc_array' plugin")

# Set up the library.
detray_add_library( detray_vc_array vc_array
   "include/detray/plugins/algebra/vc_array_definitions.hpp" )
target_link_libraries( detray_vc_array
   INTERFACE algebra::vc_vc vecmem::core )
