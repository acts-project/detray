# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::eigen' plugin")

# Set up the library.
detray_add_library( detray_eigen eigen
   "include/detray/plugins/algebra/eigen_definitions.hpp" )
target_link_libraries( detray_eigen
   INTERFACE algebra::eigen_eigen vecmem::core )
