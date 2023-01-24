# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include( detray-functions )

# Basic flags for all build modes.
foreach( mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG )
   detray_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wall" )
   detray_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wextra" )
   detray_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wno-unknown-cuda-version" )
   detray_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wshadow" )
   detray_add_flag( CMAKE_SYCL_FLAGS_${mode} "-Wunused-local-typedefs" )
endforeach()

# More rigorous tests for the Debug builds.
detray_add_flag( CMAKE_SYCL_FLAGS_DEBUG "-Werror" )
if( NOT WIN32 )
   detray_add_flag( CMAKE_SYCL_FLAGS_DEBUG "-pedantic" )
endif()

# Avoid issues coming from MSVC<->DPC++ argument differences.
if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
   foreach( mode RELEASE RELWITHDEBINFO MINSIZEREL DEBUG )
      detray_add_flag( CMAKE_SYCL_FLAGS_${mode}
         "-Wno-unused-command-line-argument" )
   endforeach()
endif()