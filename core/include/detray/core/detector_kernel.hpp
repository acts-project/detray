/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"

namespace detray {

namespace detail {

/// A functor to update the mask index in surface objects
struct mask_index_update {
    using output_type = bool;

    template <typename group_t, typename index_t, typename surface_t>
    DETRAY_HOST inline output_type operator()(const group_t& group,
                                              const index_t& /*index*/,
                                              surface_t& sf) const {
        sf.update_mask(group.size());
        return true;
    }
};

/// A functor to update the material index in surface objects
struct material_index_update {
    using output_type = bool;

    template <typename group_t, typename index_t, typename surface_t>
    DETRAY_HOST inline output_type operator()(const group_t& group,
                                              const index_t& /*index*/,
                                              surface_t& sf) const {
        sf.update_material(group.size());
        return true;
    }
};

/// A functor to perform global to local transformation
template <typename algebra_t>
struct global_to_local {

    using output_type = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        const algebra_t& trf3, const point3& pos, const vector3& dir) const {

        return mask_group[index].to_local_frame(trf3, pos, dir);
    }
};

}  // namespace detail

}  // namespace detray
