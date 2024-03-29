# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# A sanity check.
if( NOT ALGEBRA_PLUGINS_INCLUDE_VC )
   message( WARNING "Vc not available from Algebra Plugins. "
      "The configuration will likely fail." )
endif()

# Set up the library.
detray_add_library( detray_algebra_vc algebra_vc
   "include/detray/plugins/algebra/vc_array_definitions.hpp" )
target_link_libraries( detray_algebra_vc
   INTERFACE algebra::vc_cmath vecmem::core )
target_compile_definitions( detray_algebra_vc
   INTERFACE DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
             DETRAY_ALGEBRA_VC )

# Set up tests for the public header(s) of detray::algebra_vc.
detray_test_public_headers( detray_algebra_vc
   "detray/plugins/algebra/vc_array_definitions.hpp" )
