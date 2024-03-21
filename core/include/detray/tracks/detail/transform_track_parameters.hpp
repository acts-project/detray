/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/tracks/detail/track_helper.hpp"
#include "detray/tracks/free_track_parameters.hpp"

namespace detray::detail {

/// This method transforms from a bound position to a point in the global
/// 3D cartesian frame
template <typename transform3_t, typename mask_t>
DETRAY_HOST_DEVICE inline auto bound_to_free_position(
    const transform3_t& trf, const mask_t& mask,
    const typename transform3_t::point2& p,
    const typename transform3_t::vector3& dir) {

    using local_frame_t =
        typename mask_t::shape::template local_frame_type<transform3_t>;

    return local_frame_t::local_to_global(trf, mask, p, dir);
}

/// This method transforms from a global 3D cartesian position to a 2D bound
/// position in the given local coordinate frame
template <typename local_frame_t, typename transform3_t>
DETRAY_HOST_DEVICE inline auto free_to_bound_position(
    const transform3_t& trf, const typename transform3_t::point3& p,
    const typename transform3_t::vector3& dir) {

    static_assert(std::is_same_v<typename local_frame_t::loc_point,
                                 typename transform3_t::point2>,
                  "Cannot define a bound position on this shape");

    return local_frame_t::global_to_local(trf, p, dir);
}

template <typename local_frame_t, typename transform3_t>
DETRAY_HOST_DEVICE inline bound_vector<transform3_t> free_to_bound_vector(
    const transform3_t& trf3, const free_vector<transform3_t>& free_vec) {

    // Matrix operator
    using matrix_operator = typename transform3_t::matrix_actor;
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    const auto pos = track_helper().pos(free_vec);
    const auto dir = track_helper().dir(free_vec);

    const auto bound_local =
        free_to_bound_position<local_frame_t>(trf3, pos, dir);

    bound_vector<transform3_t> bound_vec;
    matrix_operator().element(bound_vec, e_bound_loc0, 0u) = bound_local[0];
    matrix_operator().element(bound_vec, e_bound_loc1, 0u) = bound_local[1];
    // The angles are defined in the global frame!
    matrix_operator().element(bound_vec, e_bound_phi, 0u) = getter::phi(dir);
    matrix_operator().element(bound_vec, e_bound_theta, 0u) =
        getter::theta(dir);
    matrix_operator().element(bound_vec, e_bound_time, 0u) =
        matrix_operator().element(free_vec, e_free_time, 0u);
    matrix_operator().element(bound_vec, e_bound_qoverp, 0u) =
        matrix_operator().element(free_vec, e_free_qoverp, 0u);

    return bound_vec;
}

template <typename transform3_t, typename mask_t>
DETRAY_HOST_DEVICE inline free_vector<transform3_t> bound_to_free_vector(
    const transform3_t& trf3, const mask_t& mask,
    const bound_vector<transform3_t>& bound_vec) {

    // Matrix operator
    using matrix_operator = typename transform3_t::matrix_actor;
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    const auto bound_local = track_helper().bound_local(bound_vec);

    const auto dir = track_helper().dir(bound_vec);

    const auto pos = bound_to_free_position(trf3, mask, bound_local, dir);

    free_vector<transform3_t> free_vec;
    matrix_operator().element(free_vec, e_free_pos0, 0u) = pos[0];
    matrix_operator().element(free_vec, e_free_pos1, 0u) = pos[1];
    matrix_operator().element(free_vec, e_free_pos2, 0u) = pos[2];
    matrix_operator().element(free_vec, e_free_time, 0u) =
        matrix_operator().element(bound_vec, e_bound_time, 0u);
    matrix_operator().element(free_vec, e_free_dir0, 0u) = dir[0];
    matrix_operator().element(free_vec, e_free_dir1, 0u) = dir[1];
    matrix_operator().element(free_vec, e_free_dir2, 0u) = dir[2];
    matrix_operator().element(free_vec, e_free_qoverp, 0u) =
        matrix_operator().element(bound_vec, e_bound_qoverp, 0u);

    return free_vec;
}

}  // namespace detray::detail
