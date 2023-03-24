# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the CMake minimum version.
cmake_minimum_required( VERSION 3.13 )

# Set up the library.
detray_add_library( detray_io_json io_json
"include/detray/io/json_defs.hpp"
"include/detray/io/json_io.hpp"
"include/detray/io/json_algebra_io.hpp"
"include/detray/io/json_geometry_io.hpp")
target_link_libraries( detray_io_json INTERFACE nlohmann_json::nlohmann_json vecmem::core )

# Make the "main I/O library" link against it.
target_link_libraries( detray_io INTERFACE detray::io_json )
