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
message( STATUS "Building Thrust as part of the Detray project" )

# Declare where to get Thrust from.
set( DETRAY_THRUST_SOURCE
   "GIT_REPOSITORY;https://github.com/NVIDIA/thrust.git;GIT_TAG;2.1.0"
   CACHE STRING "Source for Thrust, when built as part of this project" )
mark_as_advanced( DETRAY_THRUST_SOURCE )
FetchContent_Declare( Thrust ${DETRAY_THRUST_SOURCE} )

# Options used in the build of Thrust.
set( THRUST_ENABLE_INSTALL_RULES TRUE CACHE BOOL
   "Install Thrust together with this project" )

# Get it into the current directory.
FetchContent_MakeAvailable( Thrust )

# Treat the Thrust headers as "system headers", to avoid getting warnings from
# them. But not for CUDA, as that would make it pick up CUB from the CUDA
# installation. Breaking the build.
get_target_property( _incDirs _Thrust_Thrust INTERFACE_INCLUDE_DIRECTORIES )
target_include_directories( _Thrust_Thrust SYSTEM
   INTERFACE $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:${_incDirs}> )
unset( _incDirs )
