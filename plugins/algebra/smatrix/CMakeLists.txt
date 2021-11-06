# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::smatrix' plugin")

# A sanity check.
if( NOT ALGEBRA_PLUGINS_INCLUDE_SMATRIX )
   message( WARNING "SMatrix not available from Algebra Plugins. "
      "The configuration will likely fail." )
endif()

# Set up the library.
detray_add_library( detray_smatrix smatrix
   "include/detray/plugins/algebra/smatrix_definitions.hpp" )
target_link_libraries( detray_smatrix
   INTERFACE algebra::smatrix_smatrix vecmem::core )
