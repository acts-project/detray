# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the test(s) that are algebra agnostic.
detray_add_test( cpu
   "core/containers.cpp"
   "core/typed_index.cpp"
   "geometry/barcode.cpp"
   "grids/populator.cpp"
   "propagator/actor_chain.cpp"
   "utils/hash_tree.cpp"
   "utils/invalid_values.cpp"
   "utils/ranges.cpp"
   "utils/tuple_helpers.cpp"
   "utils/type_list.cpp"
   "utils/sort.cpp"
   LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::test detray::core )

# Macro setting up the CPU tests for a specific algebra plugin.
macro( detray_add_cpu_test algebra )
   # Build the test executable.
   detray_add_test( cpu_${algebra}
      "builders/detector_builder.cpp"
      "builders/grid_builder.cpp"
      "builders/homogeneous_material_builder.cpp"
      "builders/material_map_builder.cpp"
      "builders/volume_builder.cpp"
      "coordinates/cartesian2.cpp"
      "coordinates/cartesian3.cpp"
      "coordinates/cylindrical2.cpp"
      "coordinates/cylindrical3.cpp"
      "coordinates/line2.cpp"
      "coordinates/polar2.cpp"
      "core/detector.cpp"
      "core/mask_store.cpp"
      "core/transform_store.cpp"
      "detectors/telescope_detector.cpp"
      "detectors/toy_detector.cpp"
      "detectors/wire_chamber.cpp"
      "geometry/detector_surface.cpp"
      "geometry/detector_volume.cpp"
      "geometry/volume_graph.cpp"
      "grids/axis.cpp"
      "grids/grid2.cpp"
      "grids/serializer.cpp"
      "intersection/cuboid_intersector.cpp"
      "intersection/cylinder_intersector.cpp"
      "intersection/helix_intersector.cpp"
      "intersection/helix_trajectory.cpp"
      "intersection/intersection_kernel.cpp"
      "intersection/intersection2D.cpp"
      "intersection/line_intersector.cpp"
      "intersection/plane_intersector.cpp"
      "masks/annulus2D.cpp"
      "masks/cylinder.cpp"
      "masks/line.cpp"
      "masks/rectangle2D.cpp"
      "masks/ring2D.cpp"
      "masks/single3D.cpp"
      "masks/trapezoid2D.cpp"
      "masks/unbounded.cpp"
      "masks/unmasked.cpp"
      "material/energy_loss.cpp"
      "material/material_interaction.cpp"
      "material/material_maps.cpp"
      "material/materials.cpp"
      "material/stopping_power_derivative.cpp"
      "propagator/covariance_transport.cpp"
      "propagator/guided_navigator.cpp"
      "propagator/line_stepper.cpp"
      "propagator/navigator.cpp"
      "propagator/propagator.cpp"
      "propagator/rk_stepper.cpp"
      "simulation/landau_sampling.cpp"
      "simulation/particle_gun.cpp"
      "simulation/scattering.cpp"
      "simulation/track_generators.cpp"
      "surface_finders/brute_force_finder.cpp"
      "surface_finders/local_object_finder.cpp"
      "surface_finders/grid_axis.cpp"
      "surface_finders/grid_collection.cpp"
      "surface_finders/grid_populator.cpp"
      "surface_finders/grid_serializer.cpp"
      "surface_finders/grid.cpp"
      "tracks/bound_track_parameters.cpp"
      "tracks/free_track_parameters.cpp"
      "utils/bounding_volume.cpp"
      "utils/axis_rotation.cpp"
      "utils/matrix_helper.cpp"
      "utils/quadratic_equation.cpp"
      "utils/unit_vectors.cpp"
      LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core_${algebra}
                     detray::test covfie::core vecmem::core detray::io
                     detray::utils)
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
