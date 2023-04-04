# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the CMake minimum version.
cmake_minimum_required( VERSION 3.13 )

# Set up the core library.
detray_add_library( detray_io_json io_json
   "include/detray/io/json_defs.hpp"
   "include/detray/io/json_io.hpp"
   "include/detray/io/json_algebra_io.hpp"
   "include/detray/io/json_geometry_io.hpp" )
target_link_libraries( detray_io_json
   INTERFACE nlohmann_json::nlohmann_json detray::core detray::io_payload )

# Make the "main I/O library" link against it.
target_link_libraries( detray_io INTERFACE detray::io_json )

# Set up libraries using particular algebra plugins.
detray_add_library( detray_io_json_array io_json_array )
target_link_libraries( detray_io_json_array
   INTERFACE detray::io_json detray::algebra_array )

if( DETRAY_EIGEN_PLUGIN )
   detray_add_library( detray_io_json_eigen io_json_eigen )
   target_link_libraries( detray_io_json_eigen
      INTERFACE detray::io_json detray::algebra_eigen )
endif()

if( DETRAY_SMATRIX_PLUGIN )
   detray_add_library( detray_io_json_smatrix io_json_smatrix )
   target_link_libraries( detray_io_json_smatrix
      INTERFACE detray::io_json detray::algebra_smatrix )
endif()

if( DETRAY_VC_PLUGIN )
   detray_add_library( detray_io_json_vc io_json_vc )
   target_link_libraries( detray_io_json_vc
      INTERFACE detray::io_json detray::algebra_vc )
endif()

# Test the public headers with all the plugins.
set( _io_json_public_headers
   "detray/io/json_defs.hpp"
   "detray/io/json_io.hpp"
   "detray/io/json_algebra_io.hpp"
   "detray/io/json_geometry_io.hpp" )
detray_test_public_headers( detray_io_json_array
   "${_io_json_public_headers}" )
if( DETRAY_EIGEN_PLUGIN )
   detray_test_public_headers( detray_io_json_eigen
      "${_io_json_public_headers}" )
endif()
if( DETRAY_SMATRIX_PLUGIN )
   detray_test_public_headers( detray_io_json_smatrix
      "${_io_json_public_headers}" )
endif()
if( DETRAY_VC_PLUGIN )
   detray_test_public_headers( detray_io_json_vc
      "${_io_json_public_headers}" )
endif()
