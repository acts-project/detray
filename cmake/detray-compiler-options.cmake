# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include( detray-functions )

# Set the language standards to use.
set( CMAKE_CXX_STANDARD 17 CACHE STRING "The (Host) C++ standard to use" )
set( CMAKE_CUDA_STANDARD 17 CACHE STRING "The (CUDA) C++ standard to use" )

# Turn on the correct setting for the __cplusplus macro with MSVC.
if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
   detray_add_flag( CMAKE_CXX_FLAGS "/Zc:__cplusplus" )
endif()

# Turn on a number of warnings for the "known compilers".
if( ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" ) OR
    ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" ) )

   # Basic flags for all build modes.
   detray_add_flag( CMAKE_CXX_FLAGS "-Wall" )
   detray_add_flag( CMAKE_CXX_FLAGS "-Wextra" )
   detray_add_flag( CMAKE_CXX_FLAGS "-Wshadow" )
   detray_add_flag( CMAKE_CXX_FLAGS "-Wunused-local-typedefs" )

elseif( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )

   # Basic flags for all build modes.
   string( REGEX REPLACE "/W[0-9]" "" CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS}" )
   detray_add_flag( CMAKE_CXX_FLAGS "/W4" )

endif()

# Set the CUDA architecture to build code for.
set( CMAKE_CUDA_ARCHITECTURES "52" CACHE STRING
   "CUDA architectures to build device code for" )

# Make CUDA generate debug symbols for the device code as well in a debug
# build.
detray_add_flag( CMAKE_CUDA_FLAGS_DEBUG "-G" )
