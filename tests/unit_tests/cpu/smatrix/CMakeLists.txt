# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "standard" tests.
detray_add_test( smatrix
   "smatrix_actor_chain.cpp"
   "smatrix_annulus2D.cpp"
   "smatrix_axis_rotation.cpp"
   "smatrix_bounding_volume.cpp"
   "smatrix_brute_force_finder.cpp"
   "smatrix_check_simulation.cpp"
   "smatrix_container.cpp"
   "smatrix_coordinate_cartesian2.cpp"
   "smatrix_coordinate_cartesian3.cpp"
   "smatrix_coordinate_cylindrical2.cpp"
   "smatrix_coordinate_cylindrical3.cpp"
   "smatrix_coordinate_line2.cpp"
   "smatrix_coordinate_polar2.cpp"
   "smatrix_core.cpp"
   "smatrix_cuboid_intersector.cpp"
   "smatrix_cylinder_intersection.cpp"
   "smatrix_cylinder.cpp"
   "smatrix_detector.cpp"
   "smatrix_geometry_barcode.cpp"
   "smatrix_geometry_volume_graph.cpp"
   "smatrix_geometry_linking.cpp"
   "smatrix_geometry_navigation.cpp"
   "smatrix_geometry_scan.cpp"
   "smatrix_grid.cpp"
   "smatrix_grid_axis.cpp"
   "smatrix_grid_builder.cpp"
   "smatrix_grid_collection.cpp"
   "smatrix_grid_populator.cpp"
   "smatrix_grid_serializer.cpp"
   "smatrix_guided_navigator.cpp"
   "smatrix_helix_trajectory.cpp"
   "smatrix_intersection_kernel.cpp"
   "smatrix_line_intersection.cpp"
   "smatrix_line_stepper.cpp"
   "smatrix_line.cpp"
   "smatrix_mask_store.cpp"
   "smatrix_material_interaction.cpp"
   "smatrix_materials.cpp"
   "smatrix_navigator.cpp"
   "smatrix_particle_gun.cpp"
   "smatrix_path_correction.cpp"
   "smatrix_planar_intersection.cpp"
   "smatrix_propagator.cpp"
   "smatrix_rectangle2D.cpp"
   "smatrix_ring2D.cpp"
   "smatrix_single3D.cpp"
   "smatrix_stepper.cpp"
   "smatrix_telescope_detector.cpp"
   "smatrix_toy_geometry.cpp"
   "smatrix_track.cpp"
   "smatrix_track_generators.cpp"
   "smatrix_transform_store.cpp"
   "smatrix_trapezoid2D.cpp"
   "smatrix_unbounded.cpp"
   "smatrix_unmasked.cpp"
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
