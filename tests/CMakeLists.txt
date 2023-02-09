# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the common C++ flags.
include( detray-compiler-options-cpp )
include_directories( SYSTEM $<TARGET_PROPERTY:dfelibs::dfelibs,INTERFACE_INCLUDE_DIRECTORIES> )

# Include all of the code-holding sub-directories.
add_subdirectory( common )
add_subdirectory( unit_tests )
add_subdirectory( simulation )
if( DETRAY_BENCHMARKS )
   add_subdirectory( benchmarks )
endif()
