# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building Thrust as part of the Detray project" )

# Declare where to get Thrust from.
FetchContent_Declare( Thrust
   URL "https://github.com/NVIDIA/thrust/archive/refs/tags/1.15.0.tar.gz"
   URL_MD5 "01fac4fad227cfc1f455e6559dc2ab57" )

# Options used in the build of Thrust.
set( THRUST_ENABLE_HEADER_TESTING FALSE CACHE BOOL
   "Turn off Thrust's internal header tests" )
set( THRUST_ENABLE_TESTING FALSE CACHE BOOL
   "Turn off the build of the Thrust tests" )
set( THRUST_ENABLE_EXAMPLES FALSE CACHE BOOL
   "Turn off the build of the Thrust examples" )
set( THRUST_ENABLE_INSTALL_RULES TRUE CACHE BOOL
   "Install Thrust together with this project" )
set( THRUST_INSTALL_CUB_HEADERS FALSE CACHE BOOL
   "Do not try to install the CUB headers" )

# Get it into the current directory.
FetchContent_MakeAvailable( Thrust )
