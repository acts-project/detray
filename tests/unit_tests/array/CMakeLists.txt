# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "standard" tests.
detray_add_test( array
   "array_annulus2.cpp" "array_core.cpp" "array_cylinder_intersection.cpp"
   "array_cylinder3.cpp" "array_detector.cpp" "array_geometry_volume_graph.cpp"
   "array_geometry_linking.cpp" "array_geometry_navigation.cpp"
   "array_geometry_scan.cpp" "array_guided_navigator.cpp" "array_helix_gun.cpp" 
   "array_intersection_kernel.cpp" "array_line_stepper.cpp" 
   "array_mask_store.cpp" "array_navigator.cpp" "array_planar_intersection.cpp" "array_propagator.cpp" "array_rectangle2.cpp" "array_ring2.cpp" "array_rk_stepper.cpp" "array_single3.cpp" "array_surfaces_finder.cpp"
   "array_telescope_detector.cpp" "array_toy_geometry.cpp" "array_track.cpp" "array_transform_store.cpp" "array_trapezoid2.cpp"
   "array_unmasked.cpp" "array_volume.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::array )

# Set up the "conditional" tests.
if( DETRAY_IO_CSV )
   detray_add_test( array_read_detector "array_read_detector.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::io
                     detray::array )
endif()
if( DETRAY_DISPLAY )
   detray_add_test( array_display_masks "array_display_masks.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::display
                     detray::array )
endif()
