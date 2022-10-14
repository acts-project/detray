/** Detray library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/tracks/detail/track_helper.hpp"

// System include(s).
#include <climits>

namespace detray {

/** Coordinate base struct
 */
template <template <class> class Derived, typename transform3_t>
struct coordinate_base {

    /// @name Type definitions for the struct
    /// @{

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
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    /// @}

    DETRAY_HOST_DEVICE
    inline bound_vector free_to_bound_vector(
        const transform3_t& trf3, const free_vector& free_vec) const {
        const point3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        const point2 local =
            Derived<transform3_t>().global_to_local(trf3, pos, dir);

        bound_vector bound_vec;
        matrix_operator().element(bound_vec, e_bound_loc0, 0) = local[0];
        matrix_operator().element(bound_vec, e_bound_loc1, 0) = local[1];
        matrix_operator().element(bound_vec, e_bound_phi, 0) = getter::phi(dir);
        matrix_operator().element(bound_vec, e_bound_theta, 0) =
            getter::theta(dir);
        matrix_operator().element(bound_vec, e_bound_time, 0) =
            matrix_operator().element(free_vec, e_free_time, 0);
        matrix_operator().element(bound_vec, e_bound_qoverp, 0) =
            matrix_operator().element(free_vec, e_free_qoverp, 0);

        return bound_vec;
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline free_vector bound_to_free_vector(
        const transform3_t& trf3, const mask_t& mask,
        const bound_vector& bound_vec) const {

        const point2 local = track_helper().local(bound_vec);
        const vector3 dir = track_helper().dir(bound_vec);

        const auto pos =
            Derived<transform3_t>().local_to_global(trf3, mask, local, dir);

        free_vector free_vec;
        matrix_operator().element(free_vec, e_free_pos0, 0) = pos[0];
        matrix_operator().element(free_vec, e_free_pos1, 0) = pos[1];
        matrix_operator().element(free_vec, e_free_pos2, 0) = pos[2];
        matrix_operator().element(free_vec, e_free_time, 0) =
            matrix_operator().element(bound_vec, e_bound_time, 0);
        matrix_operator().element(free_vec, e_free_dir0, 0) = dir[0];
        matrix_operator().element(free_vec, e_free_dir1, 0) = dir[1];
        matrix_operator().element(free_vec, e_free_dir2, 0) = dir[2];
        matrix_operator().element(free_vec, e_free_qoverp, 0) =
            matrix_operator().element(bound_vec, e_bound_qoverp, 0);

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
        const scalar_type theta =
            matrix_operator().element(bound_vec, e_bound_theta, 0);
        const scalar_type phi =
            matrix_operator().element(bound_vec, e_bound_phi, 0);

        const scalar_type cos_theta = std::cos(theta);
        const scalar_type sin_theta = std::sin(theta);
        const scalar_type cos_phi = std::cos(phi);
        const scalar_type sin_phi = std::sin(phi);

        // Global position and direction
        const auto free_vec = bound_to_free_vector(trf3, mask, bound_vec);

        const vector3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        // Set d(x,y,z)/d(loc0, loc1)
        Derived<transform3_t>().set_bound_pos_to_free_pos_derivative(
            jac_to_global, trf3, mask, pos, dir);

        // Set d(bound time)/d(free time)
        matrix_operator().element(jac_to_global, e_free_time, e_bound_time) = 1;

        // Set d(n_x,n_y,n_z)/d(phi, theta)
        matrix_operator().element(jac_to_global, e_free_dir0, e_bound_phi) =
            -1 * sin_theta * sin_phi;
        matrix_operator().element(jac_to_global, e_free_dir0, e_bound_theta) =
            cos_theta * cos_phi;
        matrix_operator().element(jac_to_global, e_free_dir1, e_bound_phi) =
            sin_theta * cos_phi;
        matrix_operator().element(jac_to_global, e_free_dir1, e_bound_theta) =
            cos_theta * sin_phi;
        matrix_operator().element(jac_to_global, e_free_dir2, e_bound_theta) =
            -1 * sin_theta;
        matrix_operator().element(jac_to_global, e_free_qoverp,
                                  e_bound_qoverp) = 1;

        // Set d(x,y,z)/d(phi, theta)
        Derived<transform3_t>().set_bound_angle_to_free_pos_derivative(
            jac_to_global, trf3, mask, pos, dir);

        return jac_to_global;
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline free_to_bound_matrix free_to_bound_jacobian(
        const transform3_t& trf3, const mask_t& mask,
        const free_vector& free_vec) const {

        // Declare jacobian for bound to free coordinate transform
        free_to_bound_matrix jac_to_local =
            matrix_operator().template zero<e_bound_size, e_free_size>();

        // Global position and direction
        const vector3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        const scalar_type theta = getter::theta(dir);
        const scalar_type phi = getter::phi(dir);

        const scalar_type cos_theta = std::cos(theta);
        const scalar_type sin_theta = std::sin(theta);
        const scalar_type cos_phi = std::cos(phi);
        const scalar_type sin_phi = std::sin(phi);

        // Set d(loc0, loc1)/d(x,y,z)
        Derived<transform3_t>().set_free_pos_to_bound_pos_derivative(
            jac_to_local, trf3, mask, pos, dir);

        // Set d(free time)/d(bound time)
        matrix_operator().element(jac_to_local, e_bound_time, e_free_time) = 1;

        // Set d(phi, theta)/d(n_x, n_y, n_z)
        // @note This codes have a serious bug when theta is equal to zero...
        matrix_operator().element(jac_to_local, e_bound_phi, e_free_dir0) =
            -1. * sin_phi / sin_theta;
        matrix_operator().element(jac_to_local, e_bound_phi, e_free_dir1) =
            cos_phi / sin_theta;
        matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir0) =
            cos_phi * cos_theta;
        matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir1) =
            sin_phi * cos_theta;
        matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir2) =
            -1 * sin_theta;

        // Set d(Free Qop)/d(Bound Qop)
        matrix_operator().element(jac_to_local, e_bound_qoverp, e_free_qoverp) =
            1;

        return jac_to_local;
    }

    template <typename mask_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE inline free_matrix path_correction(
        const stepper_state_t& stepping, const transform3_t& trf3,
        const mask_t& mask) {

        free_matrix path_correction =
            matrix_operator().template zero<e_free_size, e_free_size>();

        // Position and direction
        const auto pos = stepping().pos();
        const auto dir = stepping().dir();

        // dir
        matrix_type<1, 3> t;
        matrix_operator().element(t, 0, 0) = dir[0];
        matrix_operator().element(t, 0, 1) = dir[1];
        matrix_operator().element(t, 0, 2) = dir[2];

        // Surface normal vector (w)
        matrix_type<1, 3> w;
        const auto normal =
            Derived<transform3_t>().normal(trf3, mask, pos, dir);
        matrix_operator().element(w, 0, 0) = normal[0];
        matrix_operator().element(w, 0, 1) = normal[1];
        matrix_operator().element(w, 0, 2) = normal[2];

        // w dot t
        const scalar_type wt = vector::dot(normal, dir);

        // transpose of t
        const matrix_type<3, 1> t_T = matrix_operator().transpose(t);

        // r correction term
        const matrix_type<1, 3> r_term = -1. / wt * w;

        // dr/dr0
        const matrix_type<3, 3> drdr0 = t_T * r_term;

        matrix_operator().template set_block<3, 3>(path_correction, drdr0,
                                                   e_free_pos0, e_free_pos0);

        if constexpr (stepper_state_t::id == stepping::id::e_rk) {
            using helix = detail::helix<transform3_t>;

            // Path length
            const auto s = stepping._s;

            const auto bvec =
                stepping._magnetic_field.at(pos[0], pos[1], pos[2]);
            const vector3 b{bvec[0], bvec[1], bvec[2]};
            // auto b = stepping._magnetic_field->get_field(pos, {});

            // helix
            helix hlx(stepping(), &b);

            // B field at the destination surface
            matrix_type<1, 3> h;
            matrix_operator().element(h, 0, 0) = hlx._h0[0];
            matrix_operator().element(h, 0, 1) = hlx._h0[1];
            matrix_operator().element(h, 0, 2) = hlx._h0[2];
            // matrix_operator().set_block(h, hlx._h0, 0, 0);

            // Normalized vector of h X t
            matrix_type<1, 3> n;
            matrix_operator().element(n, 0, 0) = hlx._n0[0];
            matrix_operator().element(n, 0, 1) = hlx._n0[1];
            matrix_operator().element(n, 0, 2) = hlx._n0[2];
            // matrix_operator().set_block(n, hlx._n0, 0, 0);

            // w cross h
            matrix_type<1, 3> wh;
            const auto _wh = vector::cross(normal, hlx._h0);
            matrix_operator().element(wh, 0, 0) = _wh[0];
            matrix_operator().element(wh, 0, 1) = _wh[1];
            matrix_operator().element(wh, 0, 2) = _wh[2];
            // matrix_operator().set_block(wh, vector::cross(normal, hlx._h0),
            // 0, 0);

            // Alpha
            const scalar_type A = hlx._alpha;

            // K
            const scalar_type K = hlx._K;

            // Ks = K*s
            const scalar_type Ks = K * s;

            // t correction term
            matrix_type<1, 3> t_term = matrix_operator().template zero<1, 3>();
            t_term = t_term +
                     (Ks - std::sin(Ks)) / K * vector::dot(hlx._h0, normal) * h;

            t_term = t_term + std::sin(Ks) / K * w;

            t_term = t_term + (1 - std::cos(Ks)) / K * wh;

            t_term = -1. / wt * t_term;

            // qoverp correction term
            const scalar_type L_term =
                -1. / wt * A * Ks / hlx.qop() * vector::dot(normal, hlx._n0);

            // dr/dt0
            const matrix_type<3, 3> drdt0 = t_T * t_term;

            // dr/dL0 (L = qoverp)
            const matrix_type<3, 1> drdL0 = t_T * L_term;

            // transpose of n
            const matrix_type<3, 1> n_T = matrix_operator().transpose(n);

            // dt/dr0
            const scalar_type AK = A * K;
            const matrix_type<3, 3> dtdr0 = AK * n_T * r_term;

            // dt/dt0
            const matrix_type<3, 3> dtdt0 = AK * n_T * t_term;

            // dt/dL0
            const matrix_type<3, 1> dtdL0 = AK * n_T * L_term;

            matrix_operator().template set_block<3, 3>(
                path_correction, drdt0, e_free_pos0, e_free_dir0);

            matrix_operator().template set_block<3, 1>(
                path_correction, drdL0, e_free_pos0, e_free_qoverp);

            matrix_operator().template set_block<3, 3>(
                path_correction, dtdr0, e_free_dir0, e_free_pos0);

            matrix_operator().template set_block<3, 3>(
                path_correction, dtdt0, e_free_dir0, e_free_dir0);

            matrix_operator().template set_block<3, 1>(
                path_correction, dtdL0, e_free_dir0, e_free_qoverp);
        }

        return path_correction;
    }
};

}  // namespace detray