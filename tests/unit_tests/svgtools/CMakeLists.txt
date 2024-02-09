# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the core tests.

include_directories( SYSTEM $<TARGET_PROPERTY:actsvg::meta,INTERFACE_INCLUDE_DIRECTORIES> )
include_directories( SYSTEM $<TARGET_PROPERTY:actsvg::core,INTERFACE_INCLUDE_DIRECTORIES> )

detray_add_test( svgtools
   "detectors.cpp"
   "grids.cpp"
   "groups.cpp"
   "intersections.cpp"
   "landmarks.cpp"
   "masks.cpp"
   "surfaces.cpp"
   "trajectories.cpp"
   "volumes.cpp"
   LINK_LIBRARIES GTest::gtest_main  detray::test detray::core_array
                  detray::utils actsvg::core actsvg::meta detray::svgtools)

detray_add_test( svgtools_web
   "web.cpp"
   LINK_LIBRARIES GTest::gtest_main detray::test detray::core_array
                  detray::utils actsvg::core actsvg::meta detray::svgtools
                  actsvg::web)
