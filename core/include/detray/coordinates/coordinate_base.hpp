/** Detray library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/tracks/detail/track_helper.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s).
#include <limits>

namespace detray {

/** Coordinate base struct
 */
template <template <class> class Derived, typename transform3_t>
struct coordinate_base {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = transform3_t;
    // Scalar type
    using scalar_type = typename transform3_t::scalar_type;
    // Point in 2D space
    using point2 = typename transform3_t::point2;
    // Point in 3D space
    using point3 = typename transform3_t::point3;
    // Vector in 3D space
    using vector3 = typename transform3_t::vector3;
    // Matrix operator
    using matrix_operator = typename transform3_t::matrix_actor;
    // Matrix size type
    using size_type = typename matrix_operator::size_ty;
    // 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    // Rotation Matrix
    using rotation_matrix = matrix_type<3, 3>;
    // Shorthand vector/matrix types related to bound track parameters.
    using bound_vector = matrix_type<e_bound_size, 1>;
    using bound_matrix = matrix_type<e_bound_size, e_bound_size>;
    // Mapping from bound track parameters.
    using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
    // Shorthand vector/matrix types related to free track parameters.
    using free_vector = matrix_type<e_free_size, 1>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;
    // Mapping from free track parameters.
    using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
    using free_to_path_matrix = matrix_type<1, e_free_size>;
    using path_to_free_matrix = matrix_type<e_free_size, 1>;
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    /// @}

    DETRAY_HOST_DEVICE
    inline bound_vector free_to_bound_vector(
        const transform3_t& trf3, const free_vector& free_vec) const {
        const point3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        const point3 local =
            Derived<transform3_t>().global_to_local(trf3, pos, dir);

        bound_vector bound_vec;
        matrix_operator().element(bound_vec, e_bound_loc0, 0u) = local[0];
        matrix_operator().element(bound_vec, e_bound_loc1, 0u) = local[1];
        // The angles are defined in the global frame!
        matrix_operator().element(bound_vec, e_bound_phi, 0u) =
            getter::phi(dir);
        matrix_operator().element(bound_vec, e_bound_theta, 0u) =
            getter::theta(dir);
        matrix_operator().element(bound_vec, e_bound_time, 0u) =
            matrix_operator().element(free_vec, e_free_time, 0u);
        matrix_operator().element(bound_vec, e_bound_qoverp, 0u) =
            matrix_operator().element(free_vec, e_free_qoverp, 0u);

        return bound_vec;
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline free_vector bound_to_free_vector(
        const transform3_t& trf3, const mask_t& mask,
        const bound_vector& bound_vec) const {

        const point2 bound_local = track_helper().bound_local(bound_vec);

        const vector3 dir = track_helper().dir(bound_vec);

        const auto pos = Derived<transform3_t>().bound_local_to_global(
            trf3, mask, bound_local, dir);

        free_vector free_vec;
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

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline bound_to_free_matrix bound_to_free_jacobian(
        const transform3_t& trf3, const mask_t& mask,
        const bound_vector& bound_vec) const {

        // Declare jacobian for bound to free coordinate transform
        bound_to_free_matrix jac_to_global =
            matrix_operator().template zero<e_free_size, e_bound_size>();

        // Get trigonometric values
        const scalar_type theta{
            matrix_operator().element(bound_vec, e_bound_theta, 0u)};
        const scalar_type phi{
            matrix_operator().element(bound_vec, e_bound_phi, 0u)};

        const scalar_type cos_theta{math::cos(theta)};
        const scalar_type sin_theta{math::sin(theta)};
        const scalar_type cos_phi{math::cos(phi)};
        const scalar_type sin_phi{math::sin(phi)};

        // Global position and direction
        const auto free_vec = bound_to_free_vector(trf3, mask, bound_vec);

        const vector3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        // Set d(x,y,z)/d(loc0, loc1)
        Derived<transform3_t>().set_bound_pos_to_free_pos_derivative(
            jac_to_global, trf3, pos, dir);

        // Set d(bound time)/d(free time)
        matrix_operator().element(jac_to_global, e_free_time, e_bound_time) =
            1.f;

        // Set d(n_x,n_y,n_z)/d(phi, theta)
        matrix_operator().element(jac_to_global, e_free_dir0, e_bound_phi) =
            -sin_theta * sin_phi;
        matrix_operator().element(jac_to_global, e_free_dir0, e_bound_theta) =
            cos_theta * cos_phi;
        matrix_operator().element(jac_to_global, e_free_dir1, e_bound_phi) =
            sin_theta * cos_phi;
        matrix_operator().element(jac_to_global, e_free_dir1, e_bound_theta) =
            cos_theta * sin_phi;
        matrix_operator().element(jac_to_global, e_free_dir2, e_bound_theta) =
            -sin_theta;
        matrix_operator().element(jac_to_global, e_free_qoverp,
                                  e_bound_qoverp) = 1.f;

        // Set d(x,y,z)/d(phi, theta)
        Derived<transform3_t>().set_bound_angle_to_free_pos_derivative(
            jac_to_global, trf3, pos, dir);

        return jac_to_global;
    }

    DETRAY_HOST_DEVICE inline free_to_bound_matrix free_to_bound_jacobian(
        const transform3_t& trf3, const free_vector& free_vec) const {

        // Declare jacobian for bound to free coordinate transform
        free_to_bound_matrix jac_to_local =
            matrix_operator().template zero<e_bound_size, e_free_size>();

        // Global position and direction
        const vector3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        const scalar_type theta{getter::theta(dir)};
        const scalar_type phi{getter::phi(dir)};

        const scalar_type cos_theta{math::cos(theta)};
        const scalar_type sin_theta{math::sin(theta)};
        const scalar_type cos_phi{math::cos(phi)};
        const scalar_type sin_phi{math::sin(phi)};

        // Set d(loc0, loc1)/d(x,y,z)
        Derived<transform3_t>().set_free_pos_to_bound_pos_derivative(
            jac_to_local, trf3, pos, dir);

        // Set d(free time)/d(bound time)
        matrix_operator().element(jac_to_local, e_bound_time, e_free_time) =
            1.f;

        // Set d(phi, theta)/d(n_x, n_y, n_z)
        // @note This codes have a serious bug when theta is equal to zero...
        matrix_operator().element(jac_to_local, e_bound_phi, e_free_dir0) =
            -sin_phi / sin_theta;
        matrix_operator().element(jac_to_local, e_bound_phi, e_free_dir1) =
            cos_phi / sin_theta;
        matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir0) =
            cos_phi * cos_theta;
        matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir1) =
            sin_phi * cos_theta;
        matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir2) =
            -sin_theta;

        // Set d(Free Qop)/d(Bound Qop)
        matrix_operator().element(jac_to_local, e_bound_qoverp, e_free_qoverp) =
            1.f;

        return jac_to_local;
    }

    DETRAY_HOST_DEVICE inline free_matrix path_correction(
        const vector3& pos, const vector3& dir, const vector3& dtds,
        const scalar dqopds, const transform3_t& trf3) const {

        free_to_path_matrix path_derivative =
            Derived<transform3_t>().path_derivative(trf3, pos, dir, dtds);

        path_to_free_matrix derivative =
            matrix_operator().template zero<e_free_size, 1u>();
        matrix_operator().element(derivative, e_free_pos0, 0u) = dir[0];
        matrix_operator().element(derivative, e_free_pos1, 0u) = dir[1];
        matrix_operator().element(derivative, e_free_pos2, 0u) = dir[2];
        matrix_operator().element(derivative, e_free_dir0, 0u) = dtds[0];
        matrix_operator().element(derivative, e_free_dir1, 0u) = dtds[1];
        matrix_operator().element(derivative, e_free_dir2, 0u) = dtds[2];
        matrix_operator().element(derivative, e_free_qoverp, 0u) = dqopds;

        return derivative * path_derivative;
    }
};

}  // namespace detray
