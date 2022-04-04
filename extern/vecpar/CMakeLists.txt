# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building vecpar as part of the Detray project" )

# Declare where to get vecpar from.
set( DETRAY_VECPAR_SOURCE
        "URL;https://github.com/wr-hamburg/vecpar/archive/refs/tags/v0.0.1.tar.gz;URL_MD5;a220b11b801a8869904fc1ad34b6aa10"
        CACHE STRING "Source for vecpar, when built as part of this project")

mark_as_advanced( DETRAY_VECPAR_SOURCE )
FetchContent_Declare( vecpar ${DETRAY_VECPAR_SOURCE} )

# Options used in the build of vecpar
set( VECPAR_SETUP_VECMEM FALSE CACHE BOOL
        "Do not set up VecMem in vecpar" )
set( VECPAR_BUILD_TESTS FALSE CACHE BOOL
        "Turn off the build of the vecpar unit tests")
set( VECPAR_BUILD_CUDA_BACKEND ${DETRAY_BUILD_CUDA} CACHE BOOL
        "Turn off the build of the vecpar unit tests")

find_package(OpenMP QUIET)
if (${OpenMP_FOUND})
    message(STATUS "Build vecpar OpenMP backend")
    set( VECPAR_BUILD_OMP_BACKEND TRUE CACHE BOOL
        "Build vecpar OpenMP backend")
endif()

if (DETRAY_BUILD_CUDA)
    find_package(CUDAToolkit QUIET)
    if (${CUDAToolkit_FOUND} AND ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
        message(STATUS "Build vecpar CUDA backend")
        set( VECPAR_BUILD_CUDA_BACKEND TRUE CACHE BOOL
            "Build vecpar CUDA backend")
    endif()
endif()

# Get it into the current directory.
FetchContent_MakeAvailable( vecpar )
