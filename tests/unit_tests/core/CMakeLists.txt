# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the core tests.
detray_add_test( core
   "grids_axis.cpp"
   "grids_grid2.cpp"
   "grids_populator.cpp"
   "grids_serializer.cpp"
   "thrust_tuple.cpp"
   "tools_hash_tree.cpp"
   "tuple_accessor.cpp"
   "utils_local_object_finder.cpp"
   "utils_ranges.cpp"
   "utils_quadratic_equation.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::core )
