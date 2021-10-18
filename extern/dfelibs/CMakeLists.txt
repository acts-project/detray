# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building dfelibs as part of the Detray project" )

# Declare where to get dfelibs from.
FetchContent_Declare( dfelibs
   URL "https://github.com/msmk0/dfelibs/archive/refs/tags/v20200416.tar.gz"
   URL_MD5 "d762f22a92d6a1aa6f88b3d48d426bad" )

# Get it into the current directory.
FetchContent_Populate( dfelibs )

# Set up a proper target and installation for the dfelibs headers. One that
# could be used correctly in the exported target(s) of Detray.
add_library( dfelibs INTERFACE )
target_include_directories( dfelibs INTERFACE
   $<BUILD_INTERFACE:${dfelibs_SOURCE_DIR}>
   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )
add_library( dfelibs::dfelibs ALIAS dfelibs )
install( DIRECTORY "${dfelibs_SOURCE_DIR}/dfe"
   DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" )
install( TARGETS dfelibs
   EXPORT dfelibs-exports )

include( CMakePackageConfigHelpers )
set( CMAKE_INSTALL_CMAKEDIR
   "${CMAKE_INSTALL_LIBDIR}/cmake/dfelibs-20200416" )
install( EXPORT dfelibs-exports
   NAMESPACE "dfelibs::"
   FILE "dfelibs-config-targets.cmake"
   DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" )
configure_package_config_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/dfelibs-config.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config.cmake"
   INSTALL_DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
   PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_CMAKEDIR
   NO_CHECK_REQUIRED_COMPONENTS_MACRO )
write_basic_package_version_file(
   "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config-version.cmake"
   VERSION "20200416"
   COMPATIBILITY "AnyNewerVersion" )
install( FILES
   "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config.cmake"
   "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/dfelibs-config-version.cmake"
   DESTINATION "${CMAKE_INSTALL_CMAKEDIR}" )
