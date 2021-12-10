# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building VecMem as part of the Detray project" )

# Declare where to get VecMem from.
set( DETRAY_VECMEM_SOURCE
   "GIT_REPOSITORY;https://github.com/acts-project/vecmem;GIT_TAG;3950939125ee08924a5b37b411e2baa846ed8682"
   CACHE STRING "Source for VecMem, when built as part of this project" )
mark_as_advanced( DETRAY_VECMEM_SOURCE )
FetchContent_Declare( VecMem ${DETRAY_VECMEM_SOURCE} )

# Options used in the build of VecMem.
set( VECMEM_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the VecMem unit tests" )

# Get it into the current directory.
FetchContent_MakeAvailable( VecMem )
