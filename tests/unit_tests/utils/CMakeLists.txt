# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the core tests.
detray_add_test( utils
   "landau_sampling.cpp"
   LINK_LIBRARIES GTest::gtest_main detray::core_array detray::utils_array )
