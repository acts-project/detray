# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the test(s) that are algebra agnostic.
detray_add_test( cpu
   "core_barcode.cpp"
   "core_container.cpp"
   "tools_actor_chain.cpp"
   LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core )

# Macro setting up the CPU tests for a specific algebra plugin.
macro( detray_add_cpu_test algebra )
   # Build the test executable.
   detray_add_test( cpu_${algebra}
      "axis_rotation.cpp"
      "check_geometry_linking.cpp"
      "check_geometry_navigation.cpp"
      "check_geometry_scan.cpp"
      "check_simulation.cpp"
      "coordinate_cartesian2.cpp"
      "coordinate_cartesian3.cpp"
      "coordinate_cylindrical2.cpp"
      "coordinate_cylindrical3.cpp"
      "coordinate_line2.cpp"
      "coordinate_polar2.cpp"
      "core_detector.cpp"
      "core_mask_store.cpp"
      "core_transform_store.cpp"
      "covariance_transport.cpp"
      "display_masks.cpp"
      "energy_loss.cpp"
      "geometry_volume.cpp"
      "geometry_volume_graph.cpp"
      "grid_axis.cpp"
      "grid_grid_builder.cpp"
      "grid_grid_collection.cpp"
      "grid_grid.cpp"
      "grid_populator.cpp"
      "grid_serializer.cpp"
      "io_read_detector.cpp"
      "line_stepper.cpp"
      "masks_annulus2D.cpp"
      "masks_cylinder.cpp"
      "masks_line.cpp"
      "masks_rectangle2D.cpp"
      "masks_ring2D.cpp"
      "masks_single3D.cpp"
      "masks_trapezoid2D.cpp"
      "masks_unbounded.cpp"
      "masks_unmasked.cpp"
      "material_interaction.cpp"
      "materials.cpp"
      "scattering.cpp"
      "sf_finder_brute_force.cpp"
      "test_core.cpp"
      "test_telescope_detector.cpp"
      "test_toy_geometry.cpp"
      "tools_bounding_volume.cpp"
      "tools_cuboid_intersector.cpp"
      "tools_cylinder_intersection.cpp"
      "tools_guided_navigator.cpp"
      "tools_helix_intersectors.cpp"
      "tools_helix_trajectory.cpp"
      "tools_intersection_kernel.cpp"
      "tools_line_intersection.cpp"
      "tools_matrix_helper.cpp"
      "tools_navigator.cpp"
      "tools_particle_gun.cpp"
      "tools_planar_intersection.cpp"
      "tools_propagator.cpp"
      "tools_stepper.cpp"
      "tools_track.cpp"
      "tools_track_generators.cpp"
      "utils_unit_vectors.cpp"
      LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core_${algebra}
                     detray::test detray_tests_common )
endmacro()

# Build the array tests.
detray_add_cpu_test( array )

# Build the Eigen tests.
if( DETRAY_EIGEN_PLUGIN )
   detray_add_cpu_test( eigen )
endif()

# Build the SMatrix tests.
if( DETRAY_SMATRIX_PLUGIN )
   detray_add_cpu_test( smatrix )
endif()

# Build the Vc tests.
if( DETRAY_VC_PLUGIN )
   detray_add_cpu_test( vc )
endif()