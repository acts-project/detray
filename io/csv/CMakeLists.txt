# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the CMake minimum version.
cmake_minimum_required( VERSION 3.13 )

# Set up the core library.
detray_add_library( detray_io_csv io_csv
  "include/detray/io/csv_io_types.hpp"
  "include/detray/io/csv_io.hpp"
  "include/detray/io/utils.hpp" )
target_link_libraries( detray_io_csv
  INTERFACE dfelibs::dfelibs vecmem::core detray::core )

# Make the "main I/O library" link against it.
target_link_libraries( detray_io INTERFACE detray::io_csv )

# Set up libraries using particular algebra plugins.
detray_add_library( detray_io_csv_array io_csv_array )
target_link_libraries( detray_io_csv_array
  INTERFACE detray::io_csv detray::algebra_array )

if( DETRAY_EIGEN_PLUGIN )
  detray_add_library( detray_io_csv_eigen io_csv_eigen )
  target_link_libraries( detray_io_csv_eigen
    INTERFACE detray::io_csv detray::algebra_eigen )
endif()

if( DETRAY_SMATRIX_PLUGIN )
  detray_add_library( detray_io_csv_smatrix io_csv_smatrix )
  target_link_libraries( detray_io_csv_smatrix
    INTERFACE detray::io_csv detray::algebra_smatrix )
endif()

if( DETRAY_VC_PLUGIN )
  detray_add_library( detray_io_csv_vc io_csv_vc )
  target_link_libraries( detray_io_csv_vc
    INTERFACE detray::io_csv detray::algebra_vc )
endif()

# Test the public headers with all the plugins.
set( _io_csv_public_headers
  "detray/io/csv_io_types.hpp"
  "detray/io/csv_io.hpp"
  "detray/io/utils.hpp" )
detray_test_public_headers( detray_io_csv_array
  "${_io_csv_public_headers}" )
if( DETRAY_EIGEN_PLUGIN )
  detray_test_public_headers( detray_io_csv_eigen
    "${_io_csv_public_headers}" )
endif()
if( DETRAY_SMATRIX_PLUGIN )
  detray_test_public_headers( detray_io_csv_smatrix
    "${_io_csv_public_headers}" )
endif()
if( DETRAY_VC_PLUGIN )
  detray_test_public_headers( detray_io_csv_vc
    "${_io_csv_public_headers}" )
endif()
