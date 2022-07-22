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

}  // namespace detray
