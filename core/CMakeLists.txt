# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::core' component")

# Set up the "core library".
file(
    GLOB _detray_core_public_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/builders/*.hpp"
    "include/detray/core/*.hpp"
    "include/detray/definitions/*.hpp"
    "include/detray/geometry/coordinates/*.hpp"
    "include/detray/geometry/shapes/*.hpp"
    "include/detray/geometry/*.hpp"
    "include/detray/grids/*.hpp"
    "include/detray/materials/*.hpp"
    "include/detray/navigation/accelerators/*.hpp"
    "include/detray/navigation/intersection/bounding_box/*.hpp"
    "include/detray/navigation/intersection/soa/*.hpp"
    "include/detray/navigation/intersection/*.hpp"
    "include/detray/navigation/*.hpp"
    "include/detray/propagator/actors/*.hpp"
    "include/detray/propagator/*.hpp"
    "include/detray/tracks/*.hpp"
    "include/detray/utils/grid/*.hpp"
    "include/detray/utils/ranges/*.hpp"
    "include/detray/utils/*.hpp"
)
file(
    GLOB _detray_core_private_headers
    "include/detray/*/detail/*.hpp"
    "include/detray/*/*/detail/*.hpp"
)
detray_add_library( detray_core core
   ${_detray_core_public_headers} ${_detray_core_private_headers}
)
target_link_libraries(detray_core INTERFACE vecmem::core)

# Generate a version header for the project.
configure_file(
    "cmake/version.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/detray/version.hpp"
    @ONLY
)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/detray/version.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/detray"
)
list(APPEND detray_core_public_headers "detray/version.hpp")
target_include_directories(
    detray_core
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}>
)

# Set up the libraries that use specific algebra plugins.
detray_add_library( detray_core_array core_array )
target_link_libraries(
    detray_core_array
    INTERFACE detray::core detray::algebra_array
)

if(DETRAY_EIGEN_PLUGIN)
    detray_add_library( detray_core_eigen core_eigen )
    target_link_libraries(
        detray_core_eigen
        INTERFACE detray::core detray::algebra_eigen
    )
endif()

if(DETRAY_SMATRIX_PLUGIN)
    detray_add_library( detray_core_smatrix core_smatrix )
    target_link_libraries(
        detray_core_smatrix
        INTERFACE detray::core detray::algebra_smatrix
    )
endif()

if(DETRAY_VC_AOS_PLUGIN)
    detray_add_library( detray_core_vc_aos core_vc_aos )
    target_link_libraries(
        detray_core_vc_aos
        INTERFACE detray::core detray::algebra_vc_aos
    )
endif()

if(DETRAY_VC_SOA_PLUGIN)
    detray_add_library( detray_core_vc_soa core_vc_soa )
    target_link_libraries(
        detray_core_vc_soa
        INTERFACE detray::core detray::algebra_vc_soa
    )
endif()

# Test the public headers of the detray core libraries.
if(BUILD_TESTING AND DETRAY_BUILD_TESTING)
    string(
        REPLACE
        "include/"
        ""
        _detray_core_public_headers
        "${_detray_core_public_headers}"
    )
    detray_test_public_headers(detray_core_array
      ${_detray_core_public_headers}
    )
    if(DETRAY_EIGEN_PLUGIN)
        detray_test_public_headers(detray_core_eigen
         ${_detray_core_public_headers}
        )
    endif()
    if(DETRAY_SMATRIX_PLUGIN)
        detray_test_public_headers(detray_core_smatrix
         ${_detray_core_public_headers}
        )
    endif()
    if(DETRAY_VC_AOS_PLUGIN)
        detray_test_public_headers(detray_core_vc_aos
         ${_detray_core_public_headers}
        )
    endif()
endif()

# Clean up.
unset(_detray_core_public_headers)
unset(_detray_core_private_headers)
