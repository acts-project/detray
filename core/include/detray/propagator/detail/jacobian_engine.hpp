/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/propagator/detail/jacobian.hpp"
#include "detray/propagator/detail/jacobian_cartesian.hpp"
#include "detray/propagator/detail/jacobian_concentric_cylindrical.hpp"
#include "detray/propagator/detail/jacobian_cylindrical.hpp"
#include "detray/propagator/detail/jacobian_line.hpp"
#include "detray/propagator/detail/jacobian_polar.hpp"
#include "detray/tracks/detail/transform_track_parameters.hpp"

namespace detray::detail {

/// @brief Generate Jacobians
template <typename algebra_t>
struct jacobian_engine {

    /// @name Type definitions for the struct
    /// @{
    using algebra_type = algebra_t;
    using transform3_type = dtransform3D<algebra_type>;
    using scalar_type = dscalar<algebra_type>;
    using point3_type = dpoint3D<algebra_type>;
    using vector3_type = dvector3D<algebra_type>;

    using bound_to_free_matrix_type = bound_to_free_matrix<algebra_type>;
    using free_to_bound_matrix_type = free_to_bound_matrix<algebra_type>;
    using free_to_path_matrix_type = free_to_path_matrix<algebra_type>;
    using path_to_free_matrix_type = path_to_free_matrix<algebra_type>;
    /// @}

    template <typename frame_t>
        requires concepts::point<typename frame_t::loc_point>
    DETRAY_HOST_DEVICE static inline void
    bound_to_free_jacobian_set_dfree_dlocal(
        bound_to_free_matrix_type& jac_to_global, const transform3_type& trf3,
        const vector3_type& pos, const vector3_type& dir) {
        using jacobian_t = jacobian<frame_t>;

        // Set d(x,y,z)/d(loc0, loc1)
        jacobian_t::set_bound_pos_to_free_pos_derivative(jac_to_global, trf3,
                                                         pos, dir);
    }

    DETRAY_HOST_DEVICE static inline void
    bound_to_free_jacobian_set_dfree_dir_dangle(
        bound_to_free_matrix_type& jac_to_global,
        const bound_parameters_vector<algebra_type>& bound_vec) {
        // Set d(bound time)/d(free time)
        getter::element(jac_to_global, e_free_time, e_bound_time) = 1.f;

        // Get trigonometric values
        const scalar_type theta{bound_vec.theta()};
        const scalar_type phi{bound_vec.phi()};
        const scalar_type cos_theta{math::cos(theta)};
        const scalar_type sin_theta{math::sin(theta)};
        const scalar_type cos_phi{math::cos(phi)};
        const scalar_type sin_phi{math::sin(phi)};

        // Set d(n_x,n_y,n_z)/d(phi, theta)
        getter::element(jac_to_global, e_free_dir0, e_bound_phi) =
            -sin_theta * sin_phi;
        getter::element(jac_to_global, e_free_dir0, e_bound_theta) =
            cos_theta * cos_phi;
        getter::element(jac_to_global, e_free_dir1, e_bound_phi) =
            sin_theta * cos_phi;
        getter::element(jac_to_global, e_free_dir1, e_bound_theta) =
            cos_theta * sin_phi;
        getter::element(jac_to_global, e_free_dir2, e_bound_theta) = -sin_theta;
        getter::element(jac_to_global, e_free_qoverp, e_bound_qoverp) = 1.f;
    }

    template <typename frame_t>
        requires concepts::point<typename frame_t::loc_point>
    DETRAY_HOST_DEVICE static inline void
    bound_to_free_jacobian_set_dfree_dangle(
        bound_to_free_matrix_type& jac_to_global, const transform3_type& trf3,
        const vector3_type& pos, const vector3_type& dir) {
        using jacobian_t = jacobian<frame_t>;

        // Set d(x,y,z)/d(phi, theta)
        jacobian_t::set_bound_angle_to_free_pos_derivative(jac_to_global, trf3,
                                                           pos, dir);
    }

    template <typename frame_t, typename mask_t>
        requires concepts::point<typename frame_t::loc_point>
    DETRAY_HOST_DEVICE static inline bound_to_free_matrix_type
    bound_to_free_jacobian(
        const transform3_type& trf3, const mask_t& mask,
        const bound_parameters_vector<algebra_type>& bound_vec) {
        bound_to_free_matrix_type jac_to_global =
            matrix::zero<bound_to_free_matrix_type>();

        // Global position and direction
        const free_track_parameters<algebra_type> free_params =
            bound_to_free_vector(trf3, mask, bound_vec);

        const vector3_type pos = free_params.pos();
        const vector3_type dir = free_params.dir();

        bound_to_free_jacobian_set_dfree_dlocal<frame_t>(jac_to_global, trf3,
                                                         pos, dir);
        bound_to_free_jacobian_set_dfree_dir_dangle(jac_to_global, bound_vec);
        bound_to_free_jacobian_set_dfree_dangle<frame_t>(jac_to_global, trf3,
                                                         pos, dir);

        return jac_to_global;
    }

    template <typename frame_t>
        requires concepts::point<typename frame_t::loc_point>
    DETRAY_HOST_DEVICE static inline void
    free_to_bound_jacobian_set_dlocal_dfree(
        free_to_bound_matrix_type& jac_to_local, const transform3_type& trf3,
        const vector3_type& pos, const vector3_type& dir) {
        using jacobian_t = jacobian<frame_t>;

        // Set d(loc0, loc1)/d(x,y,z)
        jacobian_t::set_free_pos_to_bound_pos_derivative(jac_to_local, trf3,
                                                         pos, dir);
    }

    DETRAY_HOST_DEVICE static inline void
    free_to_bound_jacobian_set_dangle_dfree_dir(
        free_to_bound_matrix_type& jac_to_local, const vector3_type& dir) {

        // Set d(free time)/d(bound time)
        getter::element(jac_to_local, e_bound_time, e_free_time) = 1.f;

        const scalar_type theta{vector::theta(dir)};
        const scalar_type phi{vector::phi(dir)};

        const scalar_type cos_theta{math::cos(theta)};
        const scalar_type sin_theta{math::sin(theta)};
        const scalar_type cos_phi{math::cos(phi)};
        const scalar_type sin_phi{math::sin(phi)};

        // Set d(phi, theta)/d(n_x, n_y, n_z)
        // @note This codes have a serious bug when theta is equal to zero...
        getter::element(jac_to_local, e_bound_phi, e_free_dir0) =
            -sin_phi / sin_theta;
        getter::element(jac_to_local, e_bound_phi, e_free_dir1) =
            cos_phi / sin_theta;
        getter::element(jac_to_local, e_bound_theta, e_free_dir0) =
            cos_phi * cos_theta;
        getter::element(jac_to_local, e_bound_theta, e_free_dir1) =
            sin_phi * cos_theta;
        getter::element(jac_to_local, e_bound_theta, e_free_dir2) = -sin_theta;

        // Set d(Free Qop)/d(Bound Qop)
        getter::element(jac_to_local, e_bound_qoverp, e_free_qoverp) = 1.f;
    }

    template <typename frame_t>
        requires concepts::point<typename frame_t::loc_point>
    DETRAY_HOST_DEVICE static inline free_to_bound_matrix_type
    free_to_bound_jacobian(
        const transform3_type& trf3,
        const free_track_parameters<algebra_type>& free_params) {

        // Declare jacobian for bound to free coordinate transform
        free_to_bound_matrix_type jac_to_local =
            matrix::zero<free_to_bound_matrix_type>();

        // Global position and direction
        const vector3_type pos = free_params.pos();
        const vector3_type dir = free_params.dir();

        free_to_bound_jacobian_set_dlocal_dfree<frame_t>(jac_to_local, trf3,
                                                         pos, dir);
        free_to_bound_jacobian_set_dangle_dfree_dir(jac_to_local, dir);

        return jac_to_local;
    }

    template <typename frame_t>
        requires concepts::point<typename frame_t::loc_point>
    DETRAY_HOST_DEVICE static inline free_to_path_matrix_type
    free_to_path_derivative(const vector3_type& pos, const vector3_type& dir,
                            const vector3_type& dtds,
                            const transform3_type& trf3) {
        using jacobian_t = jacobian<frame_t>;

        return jacobian_t::path_derivative(trf3, pos, dir, dtds);
    }

    DETRAY_HOST_DEVICE static inline path_to_free_matrix_type
    path_to_free_derivative(const vector3_type& dir, const vector3_type& dtds,
                            const scalar_type dqopds) {

        path_to_free_matrix_type derivative =
            matrix::zero<path_to_free_matrix_type>();
        getter::element(derivative, e_free_pos0, 0u) = dir[0];
        getter::element(derivative, e_free_pos1, 0u) = dir[1];
        getter::element(derivative, e_free_pos2, 0u) = dir[2];
        getter::element(derivative, e_free_dir0, 0u) = dtds[0];
        getter::element(derivative, e_free_dir1, 0u) = dtds[1];
        getter::element(derivative, e_free_dir2, 0u) = dtds[2];
        getter::element(derivative, e_free_qoverp, 0u) = dqopds;

        return derivative;
    }

    template <typename frame_t>
        requires concepts::point<typename frame_t::loc_point>
    DETRAY_HOST_DEVICE static inline free_matrix<algebra_type> path_correction(
        const vector3_type& pos, const vector3_type& dir,
        const vector3_type& dtds, const scalar_type dqopds,
        const transform3_type& trf3) {

        return path_to_free_derivative(dir, dtds, dqopds) *
               free_to_path_derivative<frame_t>(pos, dir, dtds, trf3);
    }
};

}  // namespace detray::detail
