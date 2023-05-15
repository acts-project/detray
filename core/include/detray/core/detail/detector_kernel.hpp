/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracks/tracks.hpp"

namespace detray::detail {

/// A functor to retrieve a single surface from an accelerator
struct get_surface {
    template <typename collection_t, typename index_t>
    DETRAY_HOST_DEVICE inline auto operator()(const collection_t& sf_finder,
                                              const index_t& index,
                                              const dindex sf_idx) const {

        return sf_finder[index].at(sf_idx);
    }
};

/// A functor to perform global to local transformation
template <typename algebra_t>
struct global_to_local {

    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline point3 operator()(const mask_group_t& mask_group,
                                                const index_t& index,
                                                const algebra_t& trf3,
                                                const point3& pos,
                                                const vector3& dir) const {

        return mask_group[index].to_local_frame(trf3, pos, dir);
    }
};

/// A functor to perform local to global transformation
template <typename algebra_t>
struct local_to_global {

    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline point3 operator()(const mask_group_t& mask_group,
                                                const index_t& index,
                                                const algebra_t& trf3,
                                                const point3& local) const {

        return mask_group[index].to_global_frame(trf3, local);
    }
};

/// A functor to perform free to bound vector
template <typename algebra_t>
struct free_to_bound_vector {

    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline bound_vector_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        const algebra_t& trf3, const free_vector_type& free_vec) const {

        const auto& m = mask_group[index];

        return m.local_frame().free_to_bound_vector(trf3, free_vec);
    }
};

/// A functor to perform bound to free vector
template <typename algebra_t>
struct bound_to_free_vector {

    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline free_vector_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        const algebra_t& trf3, const bound_vector_type& bound_vec) const {

        const auto& m = mask_group[index];

        return m.local_frame().bound_to_free_vector(trf3, m, bound_vec);
    }
};

}  // namespace detray::detail