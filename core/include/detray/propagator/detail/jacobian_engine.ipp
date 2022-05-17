/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

template <typename scalar_t>
detray::bound_to_free_matrix
DETRAY_HOST_DEVICE
detray::detail::jacobian_engine<scalar_t>::bound_to_free_coordinate(
    const transform3& trf3, const detray::bound_vector& bound_vec) const {

    // Declare jacobian for bound to free coordinate transform
    bound_to_free_matrix jac_to_global =
        matrix_operator().template zero<e_free_size, e_bound_size>();

    // Get trigonometric values
    const scalar_t theta =
        matrix_operator().element(bound_vec, e_bound_theta, 0);
    const scalar_t phi = matrix_operator().element(bound_vec, e_bound_phi, 0);
    const scalar_t cos_theta = std::cos(theta);
    const scalar_t sin_theta = std::sin(theta);
    const scalar_t cos_phi = std::cos(phi);
    const scalar_t sin_phi = std::sin(phi);

    // Set d(x,y,z)/d(loc0,loc1)
    const matrix_type<3, 2> bound_to_free_rotation =
        matrix_operator().template block<3, 2>(trf3.matrix(), 0, 0);

    matrix_operator().template set_block(jac_to_global, bound_to_free_rotation,
                                         e_free_pos0, e_bound_loc0);

    matrix_operator().element(jac_to_global, e_free_time, e_bound_time) = 1;
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
    matrix_operator().element(jac_to_global, e_free_qoverp, e_bound_qoverp) = 1;

    return jac_to_global;
}

template <typename scalar_t>
detray::free_to_bound_matrix
DETRAY_HOST_DEVICE
detray::detail::jacobian_engine<scalar_t>::free_to_bound_coordinate(
    const transform3& trf3, const detray::free_vector& free_vec) const {

    // Declare jacobian for free to bound coordinate transform
    free_to_bound_matrix jac_to_local =
        matrix_operator().template zero<e_bound_size, e_free_size>();

    // Free direction
    const vector3 dir = vector_engine().dir(free_vec);

    // Get trigonometric values
    const scalar_t theta = getter::theta(dir);
    const scalar_t phi = getter::phi(dir);
    const scalar_t cos_theta = std::cos(theta);
    const scalar_t sin_theta = std::sin(theta);
    const scalar_t inv_sin_theta = 1. / sin_theta;
    const scalar_t cos_phi = std::cos(phi);
    const scalar_t sin_phi = std::sin(phi);

    // Set d(loc0,loc1)/d(x,y,z)
    const matrix_type<2, 3> free_to_bound_rotation =
        matrix_operator().template block<2, 3>(trf3.matrix_inverse(), 0, 0);

    matrix_operator().template set_block(jac_to_local, free_to_bound_rotation,
                                         e_bound_loc0, e_free_pos0);

    // Set d(Free time)/d(Bound time)
    matrix_operator().element(jac_to_local, e_bound_time, e_free_time) = 1;

    // Set d(phi, theta)/d(free dir)
    matrix_operator().element(jac_to_local, e_bound_phi, e_free_dir0) =
        -1. * sin_phi * inv_sin_theta;
    matrix_operator().element(jac_to_local, e_bound_phi, e_free_dir1) =
        cos_phi * inv_sin_theta;
    matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir0) =
        cos_phi * cos_theta;
    matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir1) =
        sin_phi * cos_theta;
    matrix_operator().element(jac_to_local, e_bound_theta, e_free_dir2) =
        -1 * sin_theta;

    // Set d(Free Qop)/d(Bound Qop)
    matrix_operator().element(jac_to_local, e_bound_qoverp, e_free_qoverp) = 1;

    return jac_to_local;
}

template <typename scalar_t>
detray::free_to_path_matrix
DETRAY_HOST_DEVICE
detray::detail::jacobian_engine<scalar_t>::free_to_path_correction(
    const transform3& trf3, const free_vector& free_vec) const {

    // Declare free to path correction
    free_to_path_matrix free_to_path =
        matrix_operator().template zero<1, e_free_size>();

    // Free direction
    const vector3 dir = vector_engine().dir(free_vec);

    // The measurement frame z axis
    const matrix_type<3, 1> ref_z_axis =
        matrix_operator().template block<3, 1>(trf3.matrix(), 0, 2);

    // cosine angle between momentum direction and the measurement frame z axis
    const scalar_t dz = vector::dot(ref_z_axis, dir);

    // Correction term
    const matrix_type<1, 3> correction_term =
        -1. / dz * matrix_operator().transpose(ref_z_axis);

    matrix_operator().template set_block<1, 3>(free_to_path, correction_term, 0,
                                               e_free_pos0);

    return free_to_path;
}