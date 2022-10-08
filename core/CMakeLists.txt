# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::core' component")

# Set up the library.
detray_add_library( detray_core core
   # coordinate include(s)
   "include/detray/coordinates/cartesian2.hpp"
   "include/detray/coordinates/cartesian3.hpp"
   "include/detray/coordinates/coordinate_base.hpp"
   "include/detray/coordinates/coordinates.hpp"
   "include/detray/coordinates/cylindrical2.hpp"
   "include/detray/coordinates/cylindrical3.hpp"
   "include/detray/coordinates/line2.hpp"
   "include/detray/coordinates/polar2.hpp"
   # core include(s)
   "include/detray/core/detail/tuple_array_container.hpp"
   "include/detray/core/detail/tuple_container.hpp"
   "include/detray/core/detail/tuple_vector_container.hpp"
   "include/detray/core/detector.hpp"
   "include/detray/core/detector_kernel.hpp"
   "include/detray/core/surfaces_finder.hpp"
   "include/detray/core/transform_store.hpp"
   "include/detray/core/type_registry.hpp"
   # definitions include(s)
   "include/detray/definitions/containers.hpp"
   "include/detray/definitions/cuda_definitions.hpp"
   "include/detray/definitions/indexing.hpp"
   "include/detray/definitions/pdg_particle.hpp"
   "include/detray/definitions/qualifiers.hpp"
   "include/detray/definitions/units.hpp"
   "include/detray/definitions/detail/accessor.hpp"
   "include/detray/definitions/track_parametrization.hpp"
   # field include(s)
   "include/detray/field/constant_magnetic_field.hpp"
   # geometry include(s)
   "include/detray/geometry/surface.hpp"
   "include/detray/geometry/volume_connector.hpp"
   "include/detray/geometry/volume_graph.hpp"
   "include/detray/geometry/volume.hpp"
   # grids include(s)
   "include/detray/grids/axis.hpp"
   "include/detray/grids/grid2.hpp"
   "include/detray/grids/populator.hpp"
   "include/detray/grids/serializer2.hpp"
   # intersection include(s)
   "include/detray/intersection/concentric_cylinder_intersector.hpp"
   "include/detray/intersection/cylinder_intersector.hpp"
   "include/detray/intersection/detail/trajectories.hpp"
   "include/detray/intersection/intersection_kernel.hpp"
   "include/detray/intersection/intersection.hpp"
   "include/detray/intersection/line_intersector.hpp"
   "include/detray/intersection/plane_intersector.hpp"
   # masks include(s)
   "include/detray/masks/annulus2.hpp"
   "include/detray/masks/cylinder3.hpp"
   "include/detray/masks/line.hpp"
   "include/detray/masks/masks.hpp"
   "include/detray/masks/masks_base.hpp"
   "include/detray/masks/rectangle2.hpp"
   "include/detray/masks/ring2.hpp"
   "include/detray/masks/single3.hpp"
   "include/detray/masks/trapezoid2.hpp"
   "include/detray/masks/unmasked.hpp"
   # materials include(s)
   "include/detray/materials/detail/density_effect_data.hpp"
   "include/detray/materials/detail/relativistic_quantities.hpp"
   "include/detray/materials/interaction.hpp"
   "include/detray/materials/material.hpp"
   "include/detray/materials/material_rod.hpp"   
   "include/detray/materials/material_slab.hpp"   
   "include/detray/materials/mixture.hpp" 
   "include/detray/materials/predefined_materials.hpp"
   # propagator include(s)
   "include/detray/propagator/actors/aborters.hpp"
   "include/detray/propagator/actors/parameter_resetter.hpp"
   "include/detray/propagator/actors/parameter_transporter.hpp"
   "include/detray/propagator/actors/pointwise_material_interactor.hpp"
   "include/detray/propagator/actors/random_scatterer.hpp"
   "include/detray/propagator/actor_chain.hpp"
   "include/detray/propagator/base_actor.hpp"
   "include/detray/propagator/base_stepper.hpp"
   "include/detray/propagator/constrained_step.hpp"
   "include/detray/propagator/line_stepper.hpp"
   "include/detray/propagator/navigation_policies.hpp"
   "include/detray/propagator/navigator.hpp"
   "include/detray/propagator/propagator.hpp"
   "include/detray/propagator/rk_stepper.hpp"
   # tools include(s)
   "include/detray/tools/associator.hpp"
   "include/detray/tools/bin_association.hpp"
   "include/detray/tools/generators.hpp"
   "include/detray/tools/grid_array_helper.hpp"
   "include/detray/tools/local_object_finder.hpp"
   # tracks include(s)
   "include/detray/tracks/detail/track_helper.hpp"
   "include/detray/tracks/bound_track_parameters.hpp"
   "include/detray/tracks/free_track_parameters.hpp"
   "include/detray/tracks/tracks.hpp"
   # utils include(s)
   "include/detray/utils/axis_rotation.hpp"
   "include/detray/utils/enumerate.hpp"
   "include/detray/utils/invalid_values.hpp"
   "include/detray/utils/math.hpp"
   "include/detray/utils/matrix_helper.hpp"
   "include/detray/utils/quadratic_equation.hpp"
   "include/detray/utils/ratio.hpp"
   "include/detray/utils/type_traits.hpp" )
target_link_libraries( detray_core INTERFACE detray_io vecmem::core detray::Thrust)
