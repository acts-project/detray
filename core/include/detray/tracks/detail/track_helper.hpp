/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/tracks/detail/trigonometrics.hpp"

// System include(s).
#include <cmath>

namespace detray::detail {

template <typename matrix_actor_t>
struct track_helper {

    /// Matrix actor
    using matrix_actor = matrix_actor_t;
    /// Size type
    using size_type = typename matrix_actor_t::size_ty;
    /// Scalar type
    using scalar_type = typename matrix_actor_t::scalar_type;
    /// 2D Matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename matrix_actor::template matrix_type<ROWS, COLS>;
    /// Array type
    template <size_type N>
    using array_type = typename matrix_actor::template array_type<N>;
    /// 3-element "vector" type
    using vector3 = array_type<3>;
    /// Point in 3D space
    using point3 = vector3;
    /// Point in 2D space
    using point2 = array_type<2>;
    /// Trigonometrics
    using trigonometrics = detail::trigonometrics<scalar_type>;

    /// Track vector types
    using bound_vector = matrix_type<e_bound_size, 1>;
    using free_vector = matrix_type<e_free_size, 1>;

    DETRAY_HOST_DEVICE
    inline point3 pos(const free_vector& free_vec) const {
        return {matrix_actor().element(free_vec, e_free_pos0, 0),
                matrix_actor().element(free_vec, e_free_pos1, 0),
                matrix_actor().element(free_vec, e_free_pos2, 0)};
    }

    DETRAY_HOST_DEVICE
    inline void set_pos(free_vector& free_vec, const point3& pos) {
        matrix_actor().element(free_vec, e_free_pos0, 0) = pos[0];
        matrix_actor().element(free_vec, e_free_pos1, 0) = pos[1];
        matrix_actor().element(free_vec, e_free_pos2, 0) = pos[2];
    }

    DETRAY_HOST_DEVICE
    inline vector3 dir(const free_vector& free_vec) const {
        return {matrix_actor().element(free_vec, e_free_dir0, 0),
                matrix_actor().element(free_vec, e_free_dir1, 0),
                matrix_actor().element(free_vec, e_free_dir2, 0)};
    }

    DETRAY_HOST_DEVICE
    inline void set_dir(free_vector& free_vec, const vector3& dir) {
        matrix_actor().element(free_vec, e_free_dir0, 0) = dir[0];
        matrix_actor().element(free_vec, e_free_dir1, 0) = dir[1];
        matrix_actor().element(free_vec, e_free_dir2, 0) = dir[2];
    }

    DETRAY_HOST_DEVICE
    inline point2 local(const bound_vector& bound_vec) const {
        return {matrix_actor().element(bound_vec, e_bound_loc0, 0),
                matrix_actor().element(bound_vec, e_bound_loc1, 0)};
    }

    DETRAY_HOST_DEVICE
    inline vector3 dir(const bound_vector& bound_vec) const {
        const auto& phi = matrix_actor().element(bound_vec, e_bound_phi, 0);
        const auto& theta = matrix_actor().element(bound_vec, e_bound_theta, 0);
        const auto cosTheta = std::cos(theta);
        const auto sinTheta = std::sin(theta);

        return {std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta};
    }

    DETRAY_HOST_DEVICE
    inline trigonometrics get_trigonometrics(
        const bound_vector& bound_vec) const {
        const scalar_type theta =
            matrix_actor().element(bound_vec, e_bound_theta, 0);
        const scalar_type phi =
            matrix_actor().element(bound_vec, e_bound_phi, 0);

        return {std::cos(theta), std::sin(theta), std::cos(phi), std::sin(phi)};
    }

    DETRAY_HOST_DEVICE
    inline trigonometrics get_trigonometrics(
        const free_vector& free_vec) const {

        const vector3 dir = this->dir(free_vec);

        const scalar_type theta = getter::theta(dir);
        const scalar_type phi = getter::phi(dir);

        return {std::cos(theta), std::sin(theta), std::cos(phi), std::sin(phi)};
    }
};

}  // namespace detray::detail
