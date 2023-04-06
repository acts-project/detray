# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the host/cpu benchmarks.
add_subdirectory( cpu )

# Set up all of the "device" benchmarks.
if( DETRAY_BUILD_CUDA )
   add_subdirectory( cuda )
endif()
