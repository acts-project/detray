/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracks/detail/track_helper.hpp"

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
    // Matrix actor
    using matrix_actor = typename transform3_t::matrix_actor;
    // Matrix size type
    using size_type = typename matrix_actor::size_ty;
    // 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename matrix_actor::template matrix_type<ROWS, COLS>;
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
    using track_helper = detail::track_helper<matrix_actor>;

    /// @}

    DETRAY_HOST_DEVICE
    inline bound_vector free_to_bound_vector(
        const transform3_t& trf3, const free_vector& free_vec) const {
        const point3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        const point2 local =
            Derived<transform3_t>().global_to_local(trf3, pos, dir);

        bound_vector bound_vec;
        matrix_actor().element(bound_vec, e_bound_loc0, 0) = local[0];
        matrix_actor().element(bound_vec, e_bound_loc1, 0) = local[1];
        matrix_actor().element(bound_vec, e_bound_phi, 0) = getter::phi(dir);
        matrix_actor().element(bound_vec, e_bound_theta, 0) =
            getter::theta(dir);
        matrix_actor().element(bound_vec, e_bound_time, 0) =
            matrix_actor().element(free_vec, e_free_time, 0);
        matrix_actor().element(bound_vec, e_bound_qoverp, 0) =
            matrix_actor().element(free_vec, e_free_qoverp, 0);

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
        matrix_actor().element(free_vec, e_free_pos0, 0) = pos[0];
        matrix_actor().element(free_vec, e_free_pos1, 0) = pos[1];
        matrix_actor().element(free_vec, e_free_pos2, 0) = pos[2];
        matrix_actor().element(free_vec, e_free_time, 0) =
            matrix_actor().element(bound_vec, e_bound_time, 0);
        matrix_actor().element(free_vec, e_free_dir0, 0) = dir[0];
        matrix_actor().element(free_vec, e_free_dir1, 0) = dir[1];
        matrix_actor().element(free_vec, e_free_dir2, 0) = dir[2];
        matrix_actor().element(free_vec, e_free_qoverp, 0) =
            matrix_actor().element(bound_vec, e_bound_qoverp, 0);

        return free_vec;
    }

    template <typename mask_t>
    DETRAY_HOST_DEVICE inline bound_to_free_matrix bound_to_free_jacobian(
        const transform3_t& trf3, const mask_t& mask,
        const bound_vector& bound_vec) const {

        // Declare jacobian for bound to free coordinate transform
        bound_to_free_matrix jac_to_global =
            matrix_actor().template zero<e_free_size, e_bound_size>();

        // Get trigonometric values
        const scalar_type theta =
            matrix_actor().element(bound_vec, e_bound_theta, 0);
        const scalar_type phi =
            matrix_actor().element(bound_vec, e_bound_phi, 0);

        const scalar_type cos_theta = std::cos(theta);
        const scalar_type sin_theta = std::sin(theta);
        const scalar_type cos_phi = std::cos(phi);
        const scalar_type sin_phi = std::sin(phi);

        // Global position and direction
        const vector3 pos = track_helper().pos(bound_vec);
        const vector3 dir = track_helper().dir(bound_vec);

        // Set d(x,y,z)/d(loc0, loc1)
        Derived<transform3_t>().set_bound_pos_to_free_pos_derivative(
            jac_to_global, trf3, mask, pos, dir);

        // Set d(bound time)/d(free time)
        matrix_actor().element(jac_to_global, e_free_time, e_bound_time) = 1;

        // Set d(n_x,n_y,n_z)/d(phi, theta)
        matrix_actor().element(jac_to_global, e_free_dir0, e_bound_phi) =
            -1 * sin_theta * sin_phi;
        matrix_actor().element(jac_to_global, e_free_dir0, e_bound_theta) =
            cos_theta * cos_phi;
        matrix_actor().element(jac_to_global, e_free_dir1, e_bound_phi) =
            sin_theta * cos_phi;
        matrix_actor().element(jac_to_global, e_free_dir1, e_bound_theta) =
            cos_theta * sin_phi;
        matrix_actor().element(jac_to_global, e_free_dir2, e_bound_theta) =
            -1 * sin_theta;
        matrix_actor().element(jac_to_global, e_free_qoverp, e_bound_qoverp) =
            1;

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
            matrix_actor().template zero<e_bound_size, e_free_size>();

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
        matrix_actor().element(jac_to_local, e_bound_time, e_free_time) = 1;

        // Set d(phi, theta)/d(n_x, n_y, n_z)
        matrix_actor().element(jac_to_local, e_bound_phi, e_free_dir0) =
            -1. * sin_phi / sin_theta;
        matrix_actor().element(jac_to_local, e_bound_phi, e_free_dir1) =
            cos_phi / sin_theta;
        matrix_actor().element(jac_to_local, e_bound_theta, e_free_dir0) =
            cos_phi * cos_theta;
        matrix_actor().element(jac_to_local, e_bound_theta, e_free_dir1) =
            sin_phi * cos_theta;
        matrix_actor().element(jac_to_local, e_bound_theta, e_free_dir2) =
            -1 * sin_theta;

        // Set d(Free Qop)/d(Bound Qop)
        matrix_actor().element(jac_to_local, e_bound_qoverp, e_free_qoverp) = 1;

        return jac_to_local;
    }

    /*
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline free_to_path_matrix free_to_path_correction(
        const transform3_t& trf3, const mask_t& mask,
        const free_vector& free_vec) const {

        // Declare free to path correction
        free_to_path_matrix free_to_path =
            matrix_actor().template zero<1, e_free_size>();

        // Global position and direction
        const vector3 pos = track_helper().pos(free_vec);
        const vector3 dir = track_helper().dir(free_vec);

        // The measurement frame z axis
        const auto frame =
            Derived<transform3_t>().reference_frame(trf3, mask, pos, dir);
        const matrix_type<3, 1> ref_z_axis =
            matrix_actor()().template block<3, 1>(frame, 0, 2);

        // cosine angle between momentum direction and the measurement frame z
        // axis
        const scalar_type dz = vector::dot(ref_z_axis, dir);

        // Correction term
        const matrix_type<1, 3> correction_term =
            -1. / dz * matrix_actor().transpose(ref_z_axis);

        matrix_actor().template set_block<1, 3>(free_to_path, correction_term,
                                                0, e_free_pos0);

        return free_to_path;
    }
    */

    template <typename mask_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE inline free_matrix path_correction(
        const stepper_state_t& stepping, const transform3_t& trf3,
        const mask_t& mask) {

        /*
        using field = typename stepper_state_t::magnetic_field;
        using helix = detail::helix<transform3_t>;

        free_matrix path_correction =
            matrix_actor().template zero<e_free_size, e_free_size>();

        sd.b_first =
            _magnetic_field.get_field(stepping().pos(), context_type{});
        */
        /*
        // Direction at the surface
        const vector3& t = stepping._step_data.k4;

        // B field at the surface
        field B(stepping._step_data.b_last);
        */
        // Use Helix
        // detail::helix<transform3_t> helix hlx(trf3, &B);

        /*
        // Direction at the surface
        const vector3& t = stepping._step_data.k4;

        // B field at the surface
        const vector3& h = stepping._step_data.b_last;

        // Normalized vector of h X t
        const vector3& n = vector::normalize(vector::cross(h,t));
        */

        // dr/dr0

        // dr/dt0

        // dr/dL0 (L = qoverp)

        // dt/dr0

        // dt/dt0

        // dt/dL0

        // return path_correction;
    }
};

}  // namespace detray