# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building Algebra Plugins as part of the Detray project" )

# Declare where to get Algebra Plugins from.
FetchContent_Declare( AlgebraPlugins
   GIT_REPOSITORY "https://github.com/krasznaa/algebra-plugins.git"
   GIT_TAG "fe2a310c2f3f517fe2646e526f95db80a2109444" )

# Options used in the build of Algebra Plugins.
set( ALGEBRA_PLUGIN_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the Algebra Plugins unit tests" )
set( ALGEBRA_PLUGIN_INCLUDE_EIGEN TRUE CACHE BOOL
   "Turn on the build of algebra::eigen" )
set( ALGEBRA_PLUGIN_INCLUDE_VC TRUE CACHE BOOL
   "Turn on the build of algebra::vc_array" )
set( ALGEBRA_PLUGIN_INCLUDE_VECMEM TRUE CACHE BOOL
   "Turn on the build of algebra::vecmem_array" )

set( ALGEBRA_PLUGIN_SETUP_EIGEN3 TRUE CACHE BOOL
   "Have Algebra Plugins set up Eigen3 for itself" )
set( ALGEBRA_PLUGIN_USE_SYSTEM_EIGEN3 FALSE CACHE BOOL
   "Have Algebra Plugins pick up Eigen3 from the system" )
set( ALGEBRA_PLUGIN_SETUP_VC TRUE CACHE BOOL
   "Have Algebra Plugins set up Vc for itself" )
set( ALGEBRA_PLUGIN_USE_SYSTEM_VC FALSE CACHE BOOL
   "Have Algebra Plugins build Vc itself" )
set( ALGEBRA_PLUGIN_SETUP_VECMEM FALSE CACHE BOOL
   "Do not set up VecMem in Algebra Plugins" )
set( ALGEBRA_PLUGIN_SETUP_GOOGLETEST FALSE CACHE BOOL
   "Do not set up GoogleTest in Algebra Plugins" )

# Get it into the current directory.
FetchContent_MakeAvailable( AlgebraPlugins )
