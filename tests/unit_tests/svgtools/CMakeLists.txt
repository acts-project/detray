# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Set up the core tests.
detray_add_unit_test( svgtools
   "detectors.cpp"
   "grids.cpp"
   "groups.cpp"
   "intersections.cpp"
   "landmarks.cpp"
   "masks.cpp"
   "material.cpp"
   "surfaces.cpp"
   "trajectories.cpp"
   "volumes.cpp"
   LINK_LIBRARIES GTest::gtest_main detray::core_array detray::test_utils
                  actsvg::core actsvg::meta detray::svgtools
)

detray_add_unit_test( svgtools_web
   "web.cpp"
   LINK_LIBRARIES GTest::gtest_main detray::core_array detray::test_utils
                  actsvg::core actsvg::meta detray::svgtools actsvg::web
)
