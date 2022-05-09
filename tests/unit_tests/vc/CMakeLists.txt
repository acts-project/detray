# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "standard" tests.
detray_add_test( vc_array
   "vc_array_actor_chain.cpp" "vc_array_annulus2.cpp" "vc_array_core.cpp" 
   "vc_array_covariance_engine.cpp" "vc_array_cylinder_intersection.cpp" 
   "vc_array_cylinder3.cpp" "vc_array_detector.cpp" "vc_array_geometry_volume_graph.cpp"
   "vc_array_geometry_linking.cpp" "vc_array_geometry_navigation.cpp"
   "vc_array_geometry_scan.cpp" "vc_array_guided_navigator.cpp"
   "vc_array_helix_gun.cpp" "vc_array_intersection_kernel.cpp"
   "vc_array_mask_store.cpp" "vc_array_navigator.cpp"
   "vc_array_planar_intersection.cpp" "vc_array_propagator.cpp"
   "vc_array_rectangle2.cpp" "vc_array_ring2.cpp"
   "vc_array_single3.cpp" "vc_array_stepper.cpp"
   "vc_array_surfaces_finder.cpp" "vc_array_telescope_detector.cpp"
   "vc_array_toy_geometry.cpp" "vc_array_track.cpp"
   "vc_array_transform_store.cpp" "vc_array_trapezoid2.cpp"
   "vc_array_unmasked.cpp" "vc_array_volume.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::vc_array )

# Set up the "conditional" tests.
if( DETRAY_IO_CSV )
   detray_add_test( vc_array_read_detector "vc_array_read_detector.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::io
                     detray::vc_array )
endif()
if( DETRAY_DISPLAY )
   detray_add_test( vc_array_display_masks "vc_array_display_masks.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::display
                     detray::vc_array )
endif()
