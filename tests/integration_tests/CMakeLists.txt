# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray host integration tests")

# Set up all of the "host" tests.
if( NOT DETRAY_NO_HOST_BUILD )
   add_subdirectory( cpu )
   add_subdirectory( io )
endif()

# Set up all of the "device" tests.
if( DETRAY_BUILD_CUDA OR DETRAY_BUILD_SYCL )
   add_subdirectory( device )
endif()
