# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building Matplot++ as part of the Detray project" )

# Declare where to get Matplot++ from.
FetchContent_Declare( Matplotpp
   URL "https://github.com/alandefreitas/matplotplusplus/archive/refs/tags/v1.1.0.tar.gz"
   URL_MD5 "261061a19eda6607215bf70fa9123744" )

# Options used in the build of Matplot++.
set( BUILD_INSTALLER TRUE CACHE BOOL
   "Install Matplot++ alongside this project" )
set( BUILD_WITH_PEDANTIC_WARNINGS FALSE CACHE BOOL
   "Turn off pedantic warnings for the Matplot++ build" )

# Get it into the current directory.
FetchContent_MakeAvailable( Matplotpp )
