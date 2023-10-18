# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::io' component")

# Set up the core I/O library.
file( GLOB _detray_io_public_headers
   RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "include/detray/io/common/*.hpp"
   "include/detray/io/csv/*.hpp"
   "include/detray/io/json/*.hpp" )
detray_add_library( detray_io io
   ${_detray_io_public_headers} )
target_link_libraries( detray_io INTERFACE
   dfelibs::dfelibs nlohmann_json::nlohmann_json vecmem::core covfie::core detray::core )

# Set up libraries using particular algebra plugins.
detray_add_library( detray_io_array io_array )
target_link_libraries( detray_io_array
   INTERFACE detray::io detray::algebra_array )

if( DETRAY_EIGEN_PLUGIN )
   detray_add_library( detray_io_eigen io_eigen )
   target_link_libraries( detray_io_eigen
      INTERFACE detray::io detray::algebra_eigen )
endif()

if( DETRAY_SMATRIX_PLUGIN )
   detray_add_library( detray_io_smatrix io_smatrix )
   target_link_libraries( detray_io_smatrix
      INTERFACE detray::io detray::algebra_smatrix )
endif()

if( DETRAY_VC_PLUGIN )
   detray_add_library( detray_io_vc io_vc )
   target_link_libraries( detray_io_vc
      INTERFACE detray::io detray::algebra_vc )
endif()

# Test the public headers of the detray I/O libraries.
if( BUILD_TESTING AND DETRAY_BUILD_TESTING )
   string( REPLACE "include/" "" _detray_io_public_headers
      "${_detray_io_public_headers}" )
   detray_test_public_headers( detray_io_array
      ${_detray_io_public_headers} )
   if( DETRAY_EIGEN_PLUGIN )
      detray_test_public_headers( detray_io_eigen
         ${_detray_io_public_headers} )
   endif()
   if( DETRAY_SMATRIX_PLUGIN )
      detray_test_public_headers( detray_io_smatrix
         ${_detray_io_public_headers} )
   endif()
   if( DETRAY_VC_PLUGIN )
      detray_test_public_headers( detray_io_vc
         ${_detray_io_public_headers} )
   endif()
endif()

# Clean up.
unset( _detray_io_public_headers )
