# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "standard" tests.
detray_add_test( smatrix
   "smatrix_annulus2.cpp" "smatrix_core.cpp" "smatrix_cylinder_intersection.cpp"
   "smatrix_cylinder3.cpp" "smatrix_detector.cpp" "smatrix_geometry_graph.cpp"
   "smatrix_geometry_linking.cpp" "smatrix_geometry_navigation.cpp"
   "smatrix_geometry_scan.cpp" "smatrix_intersection_kernel.cpp"
   "smatrix_line_stepper.cpp" "smatrix_mask_store.cpp" "smatrix_navigator.cpp"
   "smatrix_planar_intersection.cpp" "smatrix_propagator.cpp"
   "smatrix_rectangle2.cpp" "smatrix_ring2.cpp" "smatrix_single3.cpp"
   "smatrix_surfaces_finder.cpp" "smatrix_toy_geometry.cpp"
   "smatrix_transform_store.cpp" "smatrix_trapezoid2.cpp" "smatrix_unmasked.cpp"
   "smatrix_volume.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::smatrix )

# Set up the "conditional" tests.
if( DETRAY_IO_CSV )
   detray_add_test( smatrix_read_detector "smatrix_read_detector.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::io
                     detray::smatrix )
endif()
if( DETRAY_DISPLAY )
   detray_add_test( smatrix_display_masks "smatrix_display_masks.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::display
                     detray::smatrix )
endif()
