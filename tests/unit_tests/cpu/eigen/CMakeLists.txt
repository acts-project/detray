# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "standard" tests.
detray_add_test( eigen
   "eigen_annulus2D.cpp"
   "eigen_brute_force_finder.cpp"
   "eigen_check_simulation.cpp"
   "eigen_coordinate_cartesian2.cpp"
   "eigen_coordinate_cartesian3.cpp"
   "eigen_coordinate_cylindrical2.cpp"
   "eigen_coordinate_cylindrical3.cpp"
   "eigen_coordinate_line2.cpp"
   "eigen_coordinate_polar2.cpp"
   "eigen_core.cpp"
   "eigen_covariance_transport.cpp"
   "eigen_cuboid_intersector.cpp"
   "eigen_cylinder_intersection.cpp"
   "eigen_cylinder.cpp"
   "eigen_detector.cpp"
   "eigen_energy_loss.cpp"
   "eigen_geometry_volume_graph.cpp"
   "eigen_geometry_linking.cpp"
   "eigen_geometry_navigation.cpp"
   "eigen_geometry_scan.cpp"
   "eigen_grid.cpp"
   "eigen_grid_axis.cpp"
   "eigen_grid_builder.cpp"
   "eigen_grid_collection.cpp"
   "eigen_grid_populator.cpp"
   "eigen_grid_serializer.cpp"
   "eigen_guided_navigator.cpp"
   "eigen_helix_intersectors.cpp"
   "eigen_helix_trajectory.cpp"
   "eigen_intersection_kernel.cpp"
   "eigen_line_intersection.cpp"
   "eigen_line_stepper.cpp"
   "eigen_line.cpp"
   "eigen_mask_store.cpp"
   "eigen_material_interaction.cpp"
   "eigen_materials.cpp"
   "eigen_matrix_helper.cpp"
   "eigen_navigator.cpp"
   "eigen_particle_gun.cpp"
   "eigen_planar_intersection.cpp"
   "eigen_projection_matrix.cpp"
   "eigen_propagator.cpp"
   "eigen_rectangle2D.cpp"
   "eigen_ring2D.cpp"
   "eigen_scattering.cpp"
   "eigen_single3D.cpp"
   "eigen_stepper.cpp"
   "eigen_telescope_detector.cpp"
   "eigen_toy_geometry.cpp"
   "eigen_track.cpp"
   "eigen_track_generators.cpp"
   "eigen_transform_store.cpp"
   "eigen_trapezoid2D.cpp"
   "eigen_unbounded.cpp"
   "eigen_unit_vectors.cpp"
   "eigen_unmasked.cpp"
   "eigen_volume.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::algebra_eigen )

detray_add_test( eigen_read_detector "eigen_read_detector.cpp"
   LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::io
                  detray::algebra_eigen )

# Set up the "conditional" test(s).
if( DETRAY_DISPLAY )
   detray_add_test( eigen_display_masks "eigen_display_masks.cpp"
      LINK_LIBRARIES GTest::gtest_main detray_tests_common detray::display
                     detray::algebra_eigen )
endif()
