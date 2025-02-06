# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Add all subdirectories.
add_subdirectory(algebra)

if(DETRAY_SVG_DISPLAY)
    add_subdirectory(svgtools)
endif()
