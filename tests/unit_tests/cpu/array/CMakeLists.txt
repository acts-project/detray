# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "standard" tests.
detray_add_test( array
   "array_actor_chain.cpp"
   "array_annulus2D.cpp"
   "array_axis_rotation.cpp"
   "array_bounding_volume.cpp"
   "array_brute_force_finder.cpp"
   "array_check_simulation.cpp"
   "array_container.cpp"
   "array_coordinate_cartesian2.cpp"
   "array_coordinate_cartesian3.cpp"
   "array_coordinate_cylindrical2.cpp"
   "array_coordinate_cylindrical3.cpp"
   "array_coordinate_line2.cpp"
   "array_coordinate_polar2.cpp"
   "array_core.cpp"
   "array_cuboid_intersector.cpp"
   "array_cylinder_intersection.cpp"
   "array_cylinder.cpp"
   "array_detector.cpp"
   "array_geometry_barcode.cpp"
   "array_geometry_volume_graph.cpp"
   "array_geometry_linking.cpp"
   "array_geometry_navigation.cpp"
   "array_geometry_scan.cpp"
   "array_grid.cpp"
   "array_grid_axis.cpp"
   "array_grid_builder.cpp"
   "array_grid_collection.cpp"
   "array_grid_populator.cpp"
   "array_grid_serializer.cpp"
   "array_guided_navigator.cpp"
   "array_helix_intersectors.cpp"
   "array_helix_trajectory.cpp"
   "array_intersection_kernel.cpp"
   "array_line_intersection.cpp"
   "array_line_stepper.cpp"
   "array_line.cpp"
   "array_mask_store.cpp"
   "array_material_interaction.cpp"
   "array_materials.cpp"
   "array_matrix_helper.cpp"
   "array_navigator.cpp"
   "array_particle_gun.cpp"
   "array_path_correction.cpp"
   "array_planar_intersection.cpp"
   "array_propagator.cpp"
   "array_rectangle2D.cpp"
   "array_ring2D.cpp"
   "array_single3D.cpp"
   "array_stepper.cpp"
   "array_telescope_detector.cpp"
   "array_toy_geometry.cpp"
   "array_track.cpp"
   "array_transform_store.cpp"
   "array_trapezoid2D.cpp"
   "array_track_generators.cpp"
   "array_unbounded.cpp"
   "array_unit_vectors.cpp"
   "array_unmasked.cpp"
   "array_volume.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::algebra_array )

detray_add_test( array_read_detector "array_read_detector.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::io
                  detray::algebra_array )

# Set up the "conditional" test(s).
if( DETRAY_DISPLAY )
   detray_add_test( array_display_masks "array_display_masks.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::display
                     detray::algebra_array )
endif()
