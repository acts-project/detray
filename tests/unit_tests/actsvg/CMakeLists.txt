# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the core tests.
detray_add_test( actsvg
   "actsvg_masks.cpp"
   "surface_svg_converter.hpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::core_array actsvg::core actsvg::meta )
