# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

message(STATUS "Building 'detray::detectors' component")

# Set the common C++ flags.
include(detray-compiler-options-cpp)

# Set up a detector library, which contains the different detector types that
# are currently available by default.
file(
    GLOB _detray_detector_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/detectors/*.hpp"
)

detray_add_library( detray_detectors detectors ${_detray_detectors_headers})
target_link_libraries(
    detray_detectors
    INTERFACE covfie::core vecmem::core detray::core detray::io
)

unset(_detray_detector_headers)
