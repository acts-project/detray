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
#include "detray/tracks/detail/track_helper.hpp"
#include "detray/tracks/free_track_parameters.hpp"

namespace detray::detail {

/// Transform a free track parameter vector to a bound track parameter vector
///
/// @param trf3 transform of the surface the bound vector should be defined on
/// @param free_param the free track parameters to be transformed
///
/// @returns the bound track parameter vector
template <typename local_frame_t>
DETRAY_HOST_DEVICE inline auto free_to_bound_vector(
    const dtransform3D<typename local_frame_t::algebra_type>& trf3,
    const free_parameters_vector<typename local_frame_t::algebra_type>&
        free_vec) {

    // Matrix operator
    using algebra_t = typename local_frame_t::algebra_type;
    using matrix_operator = dmatrix_operator<algebra_t>;

    const auto pos = free_vec.pos();
    const auto dir = free_vec.dir();

    const auto bound_local = local_frame_t::global_to_local(trf3, pos, dir);

    bound_vector<algebra_t> bound_vec;
    matrix_operator().element(bound_vec, e_bound_loc0, 0u) = bound_local[0];
    matrix_operator().element(bound_vec, e_bound_loc1, 0u) = bound_local[1];
    // The angles are defined in the global frame!
    matrix_operator().element(bound_vec, e_bound_phi, 0u) = getter::phi(dir);
    matrix_operator().element(bound_vec, e_bound_theta, 0u) =
        getter::theta(dir);
    matrix_operator().element(bound_vec, e_bound_time, 0u) = free_vec.time();
    matrix_operator().element(bound_vec, e_bound_qoverp, 0u) = free_vec.qop();

    return bound_vec;
}

/// Transform a bound track parameter vector to a free track parameter vector
///
/// @param trf3 transform of the surface the bound parameters are defined on
/// @param mask the mask of the surface the bound parameters are defined on
/// @param bound_vec the bound track vector to be transformed
///
/// @returns the free track parameter vector
template <typename mask_t>
DETRAY_HOST_DEVICE inline auto bound_to_free_vector(
    const dtransform3D<typename mask_t::algebra_type>& trf3, const mask_t& mask,
    const bound_vector<typename mask_t::algebra_type>& bound_vec) {

    // Matrix operator
    using algebra_t = typename mask_t::algebra_type;
    using local_frame_t = typename mask_t::local_frame_type;
    using matrix_operator = dmatrix_operator<algebra_t>;
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    const auto bound_local = track_helper().bound_local(bound_vec);

    const auto dir = track_helper().dir(bound_vec);

    const auto pos =
        local_frame_t::local_to_global(trf3, mask, bound_local, dir);

    // The free vector constructor expects momentum and charge, so set the
    // values explicitly instead
    free_parameters_vector<algebra_t> free_vec{};

    free_vec.set_pos(pos);
    free_vec.set_time(track_helper().time(bound_vec));
    free_vec.set_dir(dir);
    free_vec.set_qop(track_helper().qop(bound_vec));

    return free_vec;
}

}  // namespace detray::detail
