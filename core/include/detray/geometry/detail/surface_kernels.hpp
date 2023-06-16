/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray::detail {

/// A functor to perform global to local transformation
template <typename algebra_t>
struct global_to_local {

    using transform3 = algebra_t;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline point3 operator()(const mask_group_t& mask_group,
                                                const index_t& index,
                                                const transform3& trf3,
                                                const point3& pos,
                                                const vector3& dir) const {

        return mask_group[index].to_local_frame(trf3, pos, dir);
    }
};

/// A functor to perform local to global transformation
template <typename algebra_t>
struct local_to_global {

    using transform3 = algebra_t;
    using point3 = typename algebra_t::point3;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline point3 operator()(const mask_group_t& mask_group,
                                                const index_t& index,
                                                const transform3& trf3,
                                                const point3& local) const {

        return mask_group[index].to_global_frame(trf3, local);
    }
};

}  // namespace detray::detail
