# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::utils' component")

# Set up the core library.
file( GLOB _detray_utils_public_headers
   RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "include/detray/detectors/*.hpp"
   "include/detray/simulation/event_generator/*.hpp"
   "include/detray/simulation/*.hpp"
   "include/detray/tracer/renderer/*.hpp"
   "include/detray/tracer/texture/*.hpp" )
detray_add_library( detray_utils utils
    ${_detray_utils_public_headers} )
target_link_libraries( detray_utils
   INTERFACE covfie::core vecmem::core detray::core detray::io )

# Set up the libraries that use specific algebra plugins.
detray_add_library( detray_utils_array utils_array )
target_link_libraries( detray_utils_array
   INTERFACE detray::utils detray::algebra_array )

if( DETRAY_EIGEN_PLUGIN )
   detray_add_library( detray_utils_eigen utils_eigen )
   target_link_libraries( detray_utils_eigen
      INTERFACE detray::utils detray::algebra_eigen )
endif()

if( DETRAY_SMATRIX_PLUGIN )
   detray_add_library( detray_utils_smatrix utils_smatrix )
   target_link_libraries( detray_utils_smatrix
      INTERFACE detray::utils detray::algebra_smatrix )
endif()

if( DETRAY_VC_PLUGIN )
   detray_add_library( detray_utils_vc utils_vc )
   target_link_libraries( detray_utils_vc
      INTERFACE detray::utils detray::algebra_vc )
endif()

# Test the public headers of the detray utils libraries.
if( BUILD_TESTING AND DETRAY_BUILD_TESTING )
   string( REPLACE "include/" "" _detray_utils_public_headers
      "${_detray_utils_public_headers}" )
   detray_test_public_headers( detray_utils_array
      ${_detray_utils_public_headers} )
   if( DETRAY_EIGEN_PLUGIN )
      detray_test_public_headers( detray_utils_eigen
         ${_detray_utils_public_headers} )
   endif()
   if( DETRAY_SMATRIX_PLUGIN )
      detray_test_public_headers( detray_utils_smatrix
         ${_detray_utils_public_headers} )
   endif()
   if( DETRAY_VC_PLUGIN )
      detray_test_public_headers( detray_utils_vc
         ${_detray_utils_public_headers} )
   endif()
endif()

# Clean up.
unset( _detray_utils_public_headers )
