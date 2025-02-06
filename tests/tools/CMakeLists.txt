# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

include(CMakeFindDependencyMacro)

find_dependency(Boost COMPONENTS program_options)

add_library(detray_tools INTERFACE)
add_library(detray::tools ALIAS detray_tools)

target_include_directories(
    detray_tools
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(
    detray_tools
    INTERFACE Boost::program_options vecmem::core detray::io detray::test_utils
)

add_subdirectory(src)
