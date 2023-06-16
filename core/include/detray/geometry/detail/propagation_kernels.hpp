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
#include "detray/tracks/tracks.hpp"

namespace detray::detail {

/// A functor to get from a free to a bound vector
template <typename algebra_t>
struct free_to_bound_vector {

    using transform3 = algebra_t;
    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;

    // Visitor to the detector mask store that is called on the mask collection
    // that contains the mask (shape) type of the surface
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline bound_vector_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        const transform3& trf3, const free_vector_type& free_vec) const {
        // Get the concrete mask instance for this surface
        const auto& m = mask_group[index];
        // Call shape-specific behaviour on the mask
        return m.local_frame().free_to_bound_vector(trf3, free_vec);
    }
};

/// A functor to get from a bound to a free vector
template <typename algebra_t>
struct bound_to_free_vector {

    using transform3 = algebra_t;
    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline free_vector_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        const transform3& trf3, const bound_vector_type& bound_vec) const {

        const auto& m = mask_group[index];

        return m.local_frame().bound_to_free_vector(trf3, m, bound_vec);
    }
};

/// A functor to get the free-to-bound Jacobian
template <typename algebra_t>
struct free_to_bound_jacobian {

    using transform3 = algebra_t;
    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline auto operator()(
        const mask_group_t& mask_group, const index_t& index,
        const transform3& trf3, const free_vector_type& free_vec) const {

        const auto& m = mask_group[index];

        return m.local_frame().free_to_bound_jacobian(trf3, free_vec);
    }
};

/// A functor to get the bound-to-free Jacobian
template <typename algebra_t>
struct bound_to_free_jacobian {

    using transform3 = algebra_t;
    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline auto operator()(
        const mask_group_t& mask_group, const index_t& index,
        const transform3& trf3, const bound_vector_type& bound_vec) const {

        const auto& m = mask_group[index];

        return m.local_frame().bound_to_free_jacobian(trf3, m, bound_vec);
    }
};

/// A functor to get the path correction
template <typename algebra_t>
struct path_correction {

    using transform3 = algebra_t;
    using vector3 = typename algebra_t::vector3;
    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;
    using matrix_operator = typename algebra_t::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline free_matrix operator()(
        const mask_group_t& mask_group, const index_t& index,
        const transform3& trf3, const vector3& pos, const vector3& dir,
        const vector3& dtds) const {

        const auto& m = mask_group[index];

        return m.local_frame().path_correction(pos, dir, dtds, trf3);
    }
};

}  // namespace detray::detail
