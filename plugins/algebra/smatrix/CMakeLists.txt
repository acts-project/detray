# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::smatrix' plugin")

# Set up the library.
detray_add_library( detray_smatrix smatrix
   "include/detray/plugins/algebra/smatrix_definitions.hpp" )
target_link_libraries( detray_smatrix
   INTERFACE algebra::smatrix_smatrix vecmem::core )
