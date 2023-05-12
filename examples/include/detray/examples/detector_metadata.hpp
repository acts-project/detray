/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/accelerator_grid.hpp"
#include "detray/surface_finders/brute_force_finder.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/vector.hpp>

/// This example defines a detray-detector type for a cylindrical
/// silicon tracker geometry with rectangle modules for the barrel,
/// trapezoids for the endcaps and a passive cylinder for the beampipe.
/// In detray, the detector module layers are wrapped in navigation volumes,
/// which in this case are cylindrical bounding volumes. The navigation volumes
/// are attached to each other by links in their boundary surfaces (called
/// portals).
/// The different types of surfaces (e.g. sensitive modules, passive material
/// surfaces, as well as portals) can be stored in geometrical acceleration
/// data structures for a fast surface neighborhood lookup during geometry
/// navigation. The volume portals are, however, always stored in a 'brute
/// force' data structure which tests all surfaces it contains.
/// In this example detector design, volumes do not contain other volumes, so
/// the volume lookup is done using a uniform grid.
namespace detray::example {

//
// Surface Primitives, as described above
//

/// Portal link type between volumes
using nav_link = std::uint_least16_t;

/// The mask types for the detector sensitive/passive surfaces
using rectangle = mask<rectangle2D<>, nav_link>;
using trapezoid = mask<trapezoid2D<>, nav_link>;
using cylinder = mask<cylinder2D<true>, nav_link>;
// Types for portals
using cylinder_portal =
    mask<cylinder2D<false, cylinder_portal_intersector>, nav_link>;
using disc_portal = mask<ring2D<>, nav_link>;

//
// Material Description
//

/// The material types to be mapped onto the surfaces: Here homogeneous material
using slab = material_slab<detray::scalar>;

//
// Acceleration Data Structures (fast surface access during navigation)
//

// Uniform surface grid definition: bin-content: std::array<dindex, 9>
template <typename grid_shape_t, typename bin_entry_t, typename container_t>
using surface_grid_t =
    grid<coordinate_axes<grid_shape_t, false, container_t>, bin_entry_t,
         simple_serializer, regular_attacher<9>>;
// Cylindrical grid for the barrel layers
template <typename bin_entry_t, typename container_t>
using cylinder_sf_grid =
    surface_grid_t<cylinder2D<>::axes<>, bin_entry_t, container_t>;
// Disc grid for the endcap layers
template <typename bin_entry_t, typename container_t>
using disc_sf_grid = surface_grid_t<ring2D<>::axes<>, bin_entry_t, container_t>;

//
// Detector
//

/// Defines a detector that contains rectangles, trapezoids, stereo annuli,
/// passive
template <typename _bfield_backend_t =
              covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                        covfie::vector::vector_d<scalar, 3>>>
struct example_metadata {
    using bfield_backend_t = _bfield_backend_t;

    /// How to index the constituent objects in a volume
    /// If they share the same index value here, they will be added into the
    /// same acceleration data structure in every respective volume
    enum geo_objects : std::size_t {
        e_sensitive = 0,  //< sensitive module surfaces in acc data structure 0
        e_passive = 0,    //< passive material surfaces in acc data structure 0
        e_portal = 0,     //< volume portals in acc data structure 0
        e_size = 1,       //< Currently no grids, so all types of surfaces are
                          //  added to the brute force search
        e_all = e_size,
    };

    /// How a volume finds its constituent objects in the detector containers
    /// In this case: One range for sensitive/passive surfaces, oportals
    using object_link_type = dmulti_index<dindex_range, geo_objects::e_size>;

    /// How to store and link transforms. The geometry context allows to resolve
    /// the conditions data for e.g. module alignment
    template <template <typename...> class vector_t = dvector>
    using transform_store =
        single_store<transform3, vector_t, geometry_context>;

    /// Assign the mask types to the mask tuple container entries. It may be a
    /// good idea to have the most common types in the first tuple entries, in
    /// order to minimize the depth of the 'unrolling' before a mask is found
    /// in the tuple
    enum class mask_ids {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_portal_ring2 = 2,
        e_portal_cylinder2 = 3,
        e_cylinder2 = 4,
    };

    /// This is the mask collections tuple (in the detector called 'mask store')
    /// the @c regular_multi_store is a vecemem-ready tuple of vectors of
    /// the detector masks.
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_store =
        regular_multi_store<mask_ids, empty_context, tuple_t, vector_t,
                            rectangle, trapezoid, cylinder_portal, disc>;

    /// Similar to the mask store, there is a material store, which
    enum class material_ids {
        e_slab = 0,
        e_none = 1,
    };

    /// How to store and link materials. The material does not make use of
    /// conditions data ( @c empty_context )
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using material_store = regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab>;

    /// Surface descriptor type used for sensitives, passives and portals
    /// It holds the indices to the surface data in the detector data stores
    /// that were defined above
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    using source_link = dindex;
    using surface_type = surface<mask_link, material_link, transform_link,
                                 nav_link, source_link>;

    /// The acceleration data structures live in another tuple that needs to
    /// indexed correctly
    enum class sf_finder_ids {
        e_brute_force = 0,    //< test all surfaces in a volume (brute force)
        e_disc_grid = 1,      //< surface grids in the endcaps
        e_cylinder_grid = 2,  //< surface grids in the barrel
        e_default = e_brute_force,
    };

    /// The tuple store that hold the acceleration data structures for all
    /// volumes. Every collection of accelerationdata structures defines its
    /// own container and view type. Does not make use of conditions data
    /// ( @c empty_context )
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using surface_finder_store =
        multi_store<sf_finder_ids, empty_context, tuple_t,
                    brute_force_collection< surface_type, container_t> 
/*, grid_collection<disc_sf_grid<surface_type, container_t>>,
    grid_collection<cylinder_sf_grid<surface_type, container_t>>
*/>;

    /// Data structure that allows to find the current detector volume from a
    /// given position. Here: Uniform grid with a 3D cylindrical shape
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<coordinate_axes<
                 cylinder3D::axes<n_axis::bounds::e_open, n_axis::irregular,
                                  n_axis::regular, n_axis::irregular>,
                 true, container_t>,
             dindex, simple_serializer, replacer>;
};

}  // namespace detray::example
