# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "standard" tests.
detray_add_test( eigen
   "eigen_actor_chain.cpp" "eigen_annulus2.cpp" "eigen_core.cpp"
   "eigen_cylinder_intersection.cpp" "eigen_cylinder3.cpp" "eigen_detector.cpp"
   "eigen_geometry_volume_graph.cpp" "eigen_geometry_linking.cpp"
   "eigen_geometry_navigation.cpp" "eigen_geometry_scan.cpp"
   "eigen_guided_navigator.cpp" "eigen_helix_gun.cpp"
   "eigen_intersection_kernel.cpp" "eigen_mask_store.cpp" "eigen_navigator.cpp"
   "eigen_planar_intersection.cpp" "eigen_propagator.cpp" "eigen_rectangle2.cpp"
   "eigen_ring2.cpp" "eigen_single3.cpp" "eigen_stepper.cpp"
   "eigen_surfaces_finder.cpp" "eigen_telescope_detector.cpp"
   "eigen_toy_geometry.cpp" "eigen_track.cpp" "eigen_transform_store.cpp"
   "eigen_trapezoid2.cpp" "eigen_unmasked.cpp" "eigen_volume.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::eigen )

# Set up the "conditional" tests.
if( DETRAY_IO_CSV )
   detray_add_test( eigen_read_detector "eigen_read_detector.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::io
                     detray::eigen )
endif()
if( DETRAY_DISPLAY )
   detray_add_test( eigen_display_masks "eigen_display_masks.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::display
                     detray::eigen )
endif()
