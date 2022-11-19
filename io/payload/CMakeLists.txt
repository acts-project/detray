# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the CMake minimum version.
cmake_minimum_required( VERSION 3.13 )

# Set up the library.
detray_add_library( detray_io_payload io_payload
"include/detray/io/io_payload.hpp")
target_link_libraries( detray_io_payload INTERFACE)

# Make the "main I/O library" link against it.
target_link_libraries( detray_io INTERFACE detray::io_payload )
