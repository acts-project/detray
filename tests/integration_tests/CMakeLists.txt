# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "host" tests.
add_subdirectory( cpu )
add_subdirectory( io )

# Set up all of the "device" tests.
if( DETRAY_BUILD_CUDA OR DETRAY_BUILD_SYCL )
   add_subdirectory( device )
endif()
