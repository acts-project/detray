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
FetchContent_Declare( VecMem
   URL "https://github.com/acts-project/vecmem/archive/refs/tags/v0.5.0.tar.gz"
   URL_MD5 "436cb6f6fecd0bd72e67b1def6883902" )

# Options used in the build of VecMem.
set( VECMEM_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the VecMem unit tests" )

# Get it into the current directory.
FetchContent_MakeAvailable( VecMem )
