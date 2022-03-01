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
   "URL;https://github.com/nlohmann/json/archive/refs/tags/v3.10.5.tar.gz;URL_MD5;5b946f7d892fa55eabec45e76a20286b"
   CACHE STRING "Source for nlohmann, when built as part of this project" )
mark_as_advanced( DETRAY_NLOHMANN_SOURCE )
FetchContent_Declare( nlohmann_json ${DETRAY_NLOHMANN_SOURCE} )

# Set Test option sto false
set( JSON_BuildTests FALSE CACHE BOOL
   "Turn the tests in nlohmann_json on/off" )

# Set Install option to True 
set( JSON_Install TRUE CACHE BOOL
   "Turn the installation of nlohmann_json on/off" )

# Get it into the current directory.
FetchContent_MakeAvailable( nlohmann_json )