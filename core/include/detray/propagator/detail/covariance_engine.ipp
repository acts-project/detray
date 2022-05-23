/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

template <typename scalar_t>
detray::bound_matrix DETRAY_HOST_DEVICE
detray::detail::covariance_engine<scalar_t>::bound_to_bound_jacobian(
    const transform3& trf3, const detray::free_vector& free_vec,
    const bound_to_free_matrix& bound_to_free_jacobian,
    const free_matrix& free_transport_jacobian,
    const free_vector& free_to_path_derivative) const {

    // Calculate the path correction term
    free_to_path_matrix path_correction_term =
        jacobian_engine().free_to_path_correction(trf3, free_vec);

    // Calculate the free-to-bound coordinate transform jacobian
    free_to_bound_matrix free_to_bound_jacobian =
        jacobian_engine().free_to_bound_coordinate(trf3, free_vec);

    // Calculate the full jacobian from the local/bound parameters at the start
    // surface to local/bound parameters at the final surface
    // @note jac(locA->locB) = jac(gloB->locB)*(1+
    // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
    const auto I =
        matrix_operator().template identity<e_free_size, e_free_size>();

    return free_to_bound_jacobian *
           (I + free_to_path_derivative * path_correction_term) *
           free_transport_jacobian * bound_to_free_jacobian;
}

template <typename scalar_t>
DETRAY_HOST_DEVICE void
detray::detail::covariance_engine<scalar_t>::bound_to_bound_covariance_update(
    const transform3& trf3, bound_matrix& bound_covariance,
    const free_vector& free_vec,
    const bound_to_free_matrix& bound_to_free_jacobian,
    const free_matrix& free_transport_jacobian,
    const free_vector& free_to_path_derivative) const {

    // Get the full jacobian
    bound_matrix full_jacobian = this->bound_to_bound_jacobian(
        trf3, free_vec, bound_to_free_jacobian, free_transport_jacobian,
        free_to_path_derivative);

    // Update the bound covariance
    bound_covariance = full_jacobian * bound_covariance *
                       matrix_operator().transpose(full_jacobian);
}

template <typename scalar_t>
DETRAY_HOST_DEVICE void
detray::detail::covariance_engine<scalar_t>::reinitialize_jacobians(
    const transform3& trf3, const bound_vector& bound_vec,
    bound_to_free_matrix& bound_to_free_jacobian,
    free_matrix& free_transport_jacobian,
    free_vector& free_to_path_derivative) const {

    // Reset jacobian coordinate transformation at the current surface
    bound_to_free_jacobian =
        jacobian_engine().bound_to_free_coordinate(trf3, bound_vec);

    // Reset jacobian transport to identity matrix
    matrix_operator().set_identity(free_transport_jacobian);

    // Reset derivative of position and direction to zero matrix
    matrix_operator().set_zero(free_to_path_derivative);
}