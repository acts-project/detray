/** Detray library, part of the ACTS project
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

    /// Scalar type
    using scalar_type = typename transform3_t::scalar_type;
    /// Point in 2D space
    using point2 = typename transform3_t::point2;
    /// Point in 3D space
    using point3 = typename transform3_t::point3;
    /// Vector in 3D space
    using vector3 = typename transform3_t::vector3;
    /// Matrix actor
    using matrix_actor = typename transform3_t::matrix_actor;
    /// Matrix size type
    using size_type = typename matrix_actor::size_ty;
    /// 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename matrix_actor::template matrix_type<ROWS, COLS>;
    /// Shorthand vector/matrix types related to bound track parameters.
    using bound_vector = matrix_type<e_bound_size, 1>;
    using bound_matrix = matrix_type<e_bound_size, e_bound_size>;
    /// Mapping from bound track parameters.
    using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
    // Shorthand vector/matrix types related to free track parameters.
    using free_vector = matrix_type<e_free_size, 1>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;
    // Mapping from free track parameters.
    using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
    using free_to_path_matrix = matrix_type<1, e_free_size>;
    // Track helper
    using track_helper = detail::track_helper<matrix_actor>;
    // Trigonometrics
    using trigonometrics = typename track_helper::trigonometrics;

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

    DETRAY_HOST_DEVICE inline bound_to_free_matrix bound_to_free_jacobian(
        const transform3_t& trf3, const bound_vector& bound_vec) {

        // Declare jacobian for bound to free coordinate transform
        bound_to_free_matrix jac_to_global =
            matrix_actor().template zero<e_free_size, e_bound_size>();

        // Get trigonometric values
        const auto t = track_helper().get_trigonometrics(bound_vec);

        // Get d(x,y,z)/d(loc0, loc1)
        const matrix_type<3, 2> bound_to_free_rotation =
            Derived<transform3_t>().bound_to_free_rotation(trf3, t);

        matrix_actor().template set_block(jac_to_global, bound_to_free_rotation,
                                          e_free_pos0, e_bound_loc0);

        matrix_actor().element(jac_to_global, e_free_time, e_bound_time) = 1;
        matrix_actor().element(jac_to_global, e_free_dir0, e_bound_phi) =
            -1 * t.sin_theta * t.sin_phi;
        matrix_actor().element(jac_to_global, e_free_dir0, e_bound_theta) =
            t.cos_theta * t.cos_phi;
        matrix_actor().element(jac_to_global, e_free_dir1, e_bound_phi) =
            t.sin_theta * t.cos_phi;
        matrix_actor().element(jac_to_global, e_free_dir1, e_bound_theta) =
            t.cos_theta * t.sin_phi;
        matrix_actor().element(jac_to_global, e_free_dir2, e_bound_theta) =
            -1 * t.sin_theta;
        matrix_actor().element(jac_to_global, e_free_qoverp, e_bound_qoverp) =
            1;

        return jac_to_global;
    }

    DETRAY_HOST_DEVICE
    inline free_to_bound_matrix free_to_bound_jacobian(
        const transform3_t& trf3, const free_vector& free_vec) {

        // Declare jacobian for bound to free coordinate transform
        free_to_bound_matrix jac_to_local =
            matrix_actor().template zero<e_bound_size, e_free_size>();

        // Free direction
        const vector3 dir = track_helper().dir(free_vec);

        // Get trigonometric values
        const auto t = track_helper().get_trigonometrics(free_vec);

        // Get d(loc0, loc1)/d(x,y,z)
        const matrix_type<2, 3> free_to_bound_rotation =
            Derived<transform3_t>().free_to_bound_rotation(trf3, t);

        matrix_actor().template set_block(jac_to_local, free_to_bound_rotation,
                                          e_bound_loc0, e_free_pos0);

        // Set d(Free time)/d(Bound time)
        matrix_actor().element(jac_to_local, e_bound_time, e_free_time) = 1;

        // Set d(phi, theta)/d(free dir)
        matrix_actor().element(jac_to_local, e_bound_phi, e_free_dir0) =
            -1. * t.sin_phi / t.sin_theta;
        matrix_actor().element(jac_to_local, e_bound_phi, e_free_dir1) =
            t.cos_phi / t.sin_theta;
        matrix_actor().element(jac_to_local, e_bound_theta, e_free_dir0) =
            t.cos_phi * t.cos_theta;
        matrix_actor().element(jac_to_local, e_bound_theta, e_free_dir1) =
            t.sin_phi * t.cos_theta;
        matrix_actor().element(jac_to_local, e_bound_theta, e_free_dir2) =
            -1 * t.sin_theta;

        // Set d(Free Qop)/d(Bound Qop)
        matrix_actor().element(jac_to_local, e_bound_qoverp, e_free_qoverp) = 1;

        return jac_to_local;
    }
};

}  // namespace detray