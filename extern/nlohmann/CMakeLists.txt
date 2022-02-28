# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building nlohmann as part of the Detray project" )

# Declare where to get Thrust from.
set( DETRAY_NLOHMANN_SOURCE
   "GIT_REPOSITORY;git@github.com:nlohmann/json.git;GIT_TAG;v3.10.0"
   CACHE STRING "Source for nlohmann, when built as part of this project" )
mark_as_advanced( DETRAY_NLOHMANN_SOURCE )
FetchContent_Declare( nlohmann ${DETRAY_NLOHMANN_SOURCE} )

# Get it into the current directory.
FetchContent_MakeAvailable( nlohmann )