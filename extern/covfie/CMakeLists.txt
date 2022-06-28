# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Fetching covfie as part of the Detray project" )

# Declare where to get covfie from.
set( DETRAY_COVFIE_SOURCE
   "URL;https://github.com/acts-project/covfie/archive/refs/tags/v0.1.0.tar.gz;URL_MD5;e5eafed7746a86bb353b6c4ee60fa373"
   CACHE STRING "Source for covfie, when built as part of this project" )
mark_as_advanced( DETRAY_COVFIE_SOURCE )
FetchContent_Declare( covfie ${DETRAY_COVFIE_SOURCE} )

# Options used for covfie.
set( COVFIE_BUILD_EXAMPLES Off CACHE BOOL "Build covfie examples")
set( COVFIE_BUILD_TESTS Off CACHE BOOL "Build covfie tests")
set( COVFIE_BUILD_BENCHMARKS Off CACHE BOOL "Build covfie benchmarks")

set( COVFIE_PLATFORM_CPU On CACHE BOOL "Enable covfie CPU platform")
set( COVFIE_PLATFORM_CUDA On CACHE BOOL "Enable covfie CUDA platform")

set( COVFIE_REQUIRE_CXX20 Off CACHE BOOL "Enable covfie C++20 requirement")

# Get it into the current directory.
FetchContent_MakeAvailable( covfie )
