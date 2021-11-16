/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/core/detector.hpp"

namespace detray {

auto create_toy_geometry2(vecmem::memory_resource& resource) {

    // geometry type
    using geometry_t = unified_index_geometry<vecmem::vector, std::array,
                                              thrust::tuple, dindex, dindex>;

    // detector type
    using detector_t =
        detector<std::array, thrust::tuple, vecmem::vector,
                 vecmem::jagged_vector, static_transform_store<vecmem::vector>,
                 geometry_t, serializer2, std::map<dindex, std::string> >;

    // sub-geometry components type
    using surface = typename detector_t::surface;
    using mask_container = typename detector_t::mask_container;
    using transform_store = typename detector_t::transform_store;
    using disc = typename detector_t::geometry::disc;
    using cylinder = typename detector_t::geometry::cylinder;
    using rectangle = typename detector_t::geometry::rectangle;
    using surface_container =
        typename detector_t::geometry::surface_filling_container;

    // Create detector
    detector_t det("toy_geometry", resource);

    // context objects
    static_transform_store<>::context ctx0;

    // parameters
    scalar detector_half_z = 500.;
    scalar beampipe_r = 27.;
    scalar first_layer_outer_r = 38.;
    scalar second_layer_inner_r = 64.;
    scalar second_layer_outer_r = 80.;

    /**
     * beampipe volume -- volume ID = 0
     */
    auto& vol0 = det.new_volume(
        {0, beampipe_r, -1 * detector_half_z, detector_half_z, -M_PI, M_PI});

    mask_container vol0_masks(resource);
    transform_store vol0_transforms(resource);

    // Add portal/surface to beampipe
    // det.add_surface<disc>(vol0.index(), {});

    /*
    det.add_surface(det.volumes().size() - 1, surface_object, transform_object,
                    mask_object);
    */
    // negative-z disc surface
    /*
    surface_t beampipe_neg_disc(tf_store.size(), {0, m_store.group<0>().size()},
                                uni_geo.n_volumes() - 1, inv_sf_finder);
                                */
    // add mask
    // m_store.group<0>().emplace_back();

    // add transform

    // positive-z disc surface
    /*
    surface_t beampipe_pos_disc(tf_store.size(), {0, m_store.group<0>().size()},
                                uni_geo.n_volumes() - 1, inv_sf_finder);
    */
    // first layer volume -- volume ID = 1
    det.new_volume({beampipe_r, first_layer_outer_r, -1 * detector_half_z,
                    detector_half_z, -M_PI, M_PI});

    // gap volume -- volume ID = 2
    det.new_volume({first_layer_outer_r, second_layer_inner_r,
                    -1 * detector_half_z, detector_half_z, -M_PI, M_PI});

    // second layer volume -- volume ID = 3
    det.new_volume({second_layer_inner_r, second_layer_outer_r,
                    -1 * detector_half_z, detector_half_z, -M_PI, M_PI});

    return det;
}

}  // namespace detray