# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the core tests.
detray_add_test( actsvg_masks
   "actsvg_masks.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::core_array actsvg::core actsvg::meta detray::actsvg)

detray_add_test( actsvg_links
   "actsvg_links.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::core_array actsvg::core actsvg::meta detray::actsvg)

detray_add_test( actsvg_portals
   "actsvg_portals.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::core_array actsvg::core actsvg::meta detray::actsvg)

detray_add_test( actsvg_volumes
   "actsvg_volumes.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::core_array actsvg::core actsvg::meta detray::actsvg)

detray_add_test( actsvg_detectors
   "actsvg_detectors.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::core_array actsvg::core actsvg::meta detray::actsvg)

