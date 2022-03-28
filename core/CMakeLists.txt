# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::core' component")

# Set up the library.
detray_add_library( detray_core core
   # core include(s)
   "include/detray/core/detector.hpp"
   "include/detray/core/mask_store.hpp"
   "include/detray/core/surfaces_finder.hpp"
   "include/detray/core/transform_store.hpp"
   "include/detray/core/type_registry.hpp"
   # definitions include(s)
   "include/detray/definitions/containers.hpp"
   "include/detray/definitions/cuda_definitions.hpp"
   "include/detray/definitions/indexing.hpp"  
   "include/detray/definitions/qualifiers.hpp"
   "include/detray/definitions/track_parameterization.hpp"
   "include/detray/definitions/units.hpp"
   "include/detray/definitions/detail/accessor.hpp"  
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
   "include/detray/intersection/intersection_kernel.hpp"
   "include/detray/intersection/intersection.hpp"
   "include/detray/intersection/planar_intersector.hpp"
   "include/detray/intersection/quadratic_equation.hpp"   
   "include/detray/intersection/unbound.hpp"   
   # masks include(s)
   "include/detray/masks/annulus2.hpp"
   "include/detray/masks/cylinder3.hpp"
   "include/detray/masks/masks.hpp"
   "include/detray/masks/rectangle2.hpp"
   "include/detray/masks/ring2.hpp"
   "include/detray/masks/single3.hpp"
   "include/detray/masks/trapezoid2.hpp"
   "include/detray/masks/unmasked.hpp"   
   # propagator include(s)
   "include/detray/propagator/base_stepper.hpp"
   "include/detray/propagator/line_stepper.hpp"
   "include/detray/propagator/navigator.hpp"
   "include/detray/propagator/propagator.hpp"
   "include/detray/propagator/rk_stepper.hpp"
   "include/detray/propagator/track.hpp"   
   # tools include(s)
   "include/detray/tools/associator.hpp"
   "include/detray/tools/bin_association.hpp"
   "include/detray/tools/generators.hpp"
   "include/detray/tools/grid_array_helper.hpp"
   "include/detray/tools/local_object_finder.hpp"
   # utils include(s)   
   "include/detray/utils/enumerate.hpp"
   "include/detray/utils/invalid_values.hpp" )
target_link_libraries( detray_core INTERFACE vecmem::core detray::Thrust)
