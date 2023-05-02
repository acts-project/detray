# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Silence FetchContent warnings with CMake >=3.24.
if( POLICY CMP0135 )
   cmake_policy( SET CMP0135 NEW )
endif()

# Tell the user what's happening.
message( STATUS "Building dfelibs as part of the Detray project" )

# Declare where to get dfelibs from.
set( DETRAY_DFELIBS_SOURCE
   "URL;https://github.com/acts-project/dfelibs/archive/refs/tags/v20211029.tar.gz;URL_MD5;87fb09c5a11b98250f5e266e9cd501ea"
   CACHE STRING "Source for dfelibs, when built as part of this project" )
mark_as_advanced( DETRAY_DFELIBS_SOURCE )
FetchContent_Declare( dfelibs ${DETRAY_DFELIBS_SOURCE} )

# Options used in the build of dfelibs.
set( dfelibs_BUILD_EXAMPLES FALSE CACHE BOOL
   "Turn off the build of the dfelibs examples" )
set( dfelibs_BUILD_UNITTESTS FALSE CACHE BOOL
   "Turn off the build of the dfelibs unit tests" )
set( dfelibs_ENABLE_INSTALL TRUE CACHE BOOL
   "Have dfelibs be installed together with the rest of the project" )

# Get it into the current directory.
FetchContent_MakeAvailable( dfelibs )
