/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/concepts/algebra.hpp"
#include "detray/concepts/index.hpp"
#include "detray/concepts/mask_group.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray::detail {

/// A functor to retrieve a single surface from an accelerator
struct get_surface {
    template <typename collection_t, CONSTRAINT(concepts::index) index_t>
    DETRAY_HOST_DEVICE inline auto operator()(const collection_t& sf_finder,
                                              const index_t& index,
                                              const dindex sf_idx) const {

        return sf_finder[index].at(sf_idx);
    }
};

/// A functor to perform global to local transformation
template <CONSTRAINT(concepts::algebra) algebra_t>
struct global_to_local {

    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    template <CONSTRAINT(concepts::mask_group) mask_group_t,
              CONSTRAINT(concepts::index) index_t>
    DETRAY_HOST_DEVICE inline auto operator()(const mask_group_t& mask_group,
                                              const index_t& index,
                                              const algebra_t& trf3,
                                              const point3& pos,
                                              const vector3& dir) const {

        return mask_group[index].to_local_frame(trf3, pos, dir);
    }
};

}  // namespace detray::detail
