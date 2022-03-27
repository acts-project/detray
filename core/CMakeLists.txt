# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::core' component")

# Set up the library.
detray_add_library( detray_core core
   "include/detray/core/detector.hpp"
   "include/detray/core/intersection.hpp"
   "include/detray/core/mask_store.hpp"
   "include/detray/core/surfaces_finder.hpp"
   "include/detray/core/transform_store.hpp"
   "include/detray/core/type_registry.hpp"
   "include/detray/definitions/detail/accessor.hpp"
   "include/detray/definitions/cuda_defs.hpp"
   "include/detray/definitions/qualifiers.hpp"
   "include/detray/definitions/invalid_values.hpp"
   "include/detray/definitions/units.hpp"
   "include/detray/field/constant_magnetic_field.hpp"
   "include/detray/geometry/index_geometry.hpp"
   "include/detray/geometry/surface.hpp"
   "include/detray/geometry/volume_connector.hpp"
   "include/detray/geometry/volume.hpp"
   "include/detray/grids/associator.hpp"
   "include/detray/grids/axis.hpp"
   "include/detray/grids/grid2.hpp"
   "include/detray/grids/populator.hpp"
   "include/detray/grids/serializer2.hpp"
   "include/detray/masks/annulus2.hpp"
   "include/detray/masks/cylinder3.hpp"
   "include/detray/masks/masks.hpp"
   "include/detray/masks/rectangle2.hpp"
   "include/detray/masks/ring2.hpp"
   "include/detray/masks/single3.hpp"
   "include/detray/masks/trapezoid2.hpp"
   "include/detray/masks/unmasked.hpp"
   "include/detray/tools/base_stepper.hpp"
   "include/detray/tools/bin_association.hpp"
   "include/detray/tools/concentric_cylinder_intersector.hpp"
   "include/detray/toolt/constrained_step.hpp"
   "include/detray/tools/cylinder_intersector.hpp"
   "include/detray/tools/volume_graph.hpp"
   "include/detray/tools/grid_array_helper.hpp"
   "include/detray/tools/intersection_kernel.hpp"
   "include/detray/tools/line_stepper.hpp"
   "include/detray/tools/local_object_finder.hpp"
   "include/detray/tools/navigator.hpp"
   "include/detray/tools/planar_intersector.hpp"
   "include/detray/tools/propagator.hpp"
   "include/detray/tools/track.hpp"   
   "include/detray/tools/track_parameterization.hpp"
   "include/detray/tools/rk_stepper.hpp"
   "include/detray/utils/containers.hpp"
   "include/detray/utils/enumerate.hpp"
   "include/detray/utils/generators.hpp"
   "include/detray/utils/indexing.hpp"
   "include/detray/utils/quadratic_equation.hpp"
   "include/detray/utils/ray_gun.hpp"
   "include/detray/utils/unbound.hpp" )
target_link_libraries( detray_core INTERFACE vecmem::core detray::Thrust)
