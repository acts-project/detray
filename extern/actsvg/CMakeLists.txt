# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Silence FetchContent warnings with CMake >=3.24.
if( POLICY CMP0135 )
   cmake_policy( SET CMP0135 NEW )
endif()

# Tell the user what's happening.
message( STATUS "Building actsvg as part of the Detray project" )

# Declare where to get Actsvg from.
set( DETRAY_ACTSVG_GIT_REPOSITORY "https://github.com/acts-project/actsvg.git"
   CACHE STRING "Git repository to take actsvg from" )
set( DETRAY_ACTSVG_GIT_TAG "v0.4.44" CACHE STRING "Version of actsvg to build" )
mark_as_advanced( DETRAY_ACTSVG_GIT_REPOSITORY DETRAY_ACTSVG_GIT_TAG )
FetchContent_Declare( actsvg
   GIT_REPOSITORY "${DETRAY_ACTSVG_GIT_REPOSITORY}"
   GIT_TAG "${DETRAY_ACTSVG_GIT_TAG}" )

# Make sure the web plugin of ACTSVG is built
set(ACTSVG_BUILD_WEB ON CACHE BOOL "Build the web plugin of ACTSVG")

# Now set up its build.
FetchContent_MakeAvailable( actsvg )

# Treat the actvg headers as "system headers", to avoid getting warnings from
# them.
get_target_property( _incDirs actsvg_core INTERFACE_INCLUDE_DIRECTORIES )
target_include_directories( actsvg_core SYSTEM INTERFACE ${_incDirs} )
unset( _incDirs )

get_target_property( _incDirs actsvg_meta INTERFACE_INCLUDE_DIRECTORIES )
target_include_directories( actsvg_meta SYSTEM INTERFACE ${_incDirs} )
unset( _incDirs )
