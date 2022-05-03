/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parameterization.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"

namespace detray {

namespace detail {

// Covariance engine agnostic to surface type
template <typename scalar_t>
struct covariance_engine {

    using jacobian_engine = detail::jacobian_engine<scalar_t>;
    using vector_engine = typename jacobian_engine::vector_engine;
    using transform3 = typename vector_engine::transform3;
    using matrix_operator = typename jacobian_engine::matrix_operator;

    /** Function to get the full jacobian for surface-to-surface propagation
     *
     * @param trf3 is the transform matrix of the destination surface
     * @param free_vec is the input free vector at the destination surface
     * @param bound_to_free_jacobian is the coordinate transform jacobian
     * at departure surface.
     * @param free_transport_jacobian is the transport jacobian
     * @param free_to_path_derivative is the derivative of free parameter w.r.t
     * path
     * @returns bound to bound jacobian
     */
    DETRAY_HOST_DEVICE inline bound_matrix bound_to_bound_jacobian(
        const transform3& trf3, const free_vector& free_vec,
        const bound_to_free_matrix& bound_to_free_jacobian,
        const free_matrix& free_transport_jacobian,
        const free_vector& free_to_path_derivative) const;

    /** Function to update the covariance matrix for surface-to-surface
     * propagation
     *
     * @param trf3 is the transform matrix of the destination surface
     * @param bound_covariance is the covariance at the departure surface, which
     * is to be updated for the destination surface
     * @param free_vec is the input free vector at the destination surface
     * @param bound_to_free_jacobian is the coordinate transform jacobian
     * at the departure surface.
     * @param free_transport_jacobian is the transport jacobian
     * @param free_to_path_derivative is the derivative of free parameter w.r.t
     * path at the destination surface
     */
    DETRAY_HOST_DEVICE inline void bound_to_bound_covariance_update(
        const transform3& trf3, bound_matrix& bound_covariance,
        const free_vector& free_vec,
        const bound_to_free_matrix& bound_to_free_jacobian,
        const free_matrix& free_transport_jacobian,
        const free_vector& free_to_path_derivative) const;

    /** Function to reinitialize the jacobian after updating covariance
     *
     * @param trf3 is the transform matrix of the destination surface
     * @param bound_vec is the input bound vector at the destination surface
     * @param bound_to_free_jacobian is the coordinate transform jacobian
     * at the destination surface
     * @param free_transport_jacobian is the transport jacobian
     * @param free_to_path_derivative is the derivative of free parameter w.r.t
     * path at the destination surface
     */
    DETRAY_HOST_DEVICE inline void reinitialize_jacobians(
        const transform3& trf3, const bound_vector& bound_vec,
        bound_to_free_matrix& bound_to_free_jacobian,
        free_matrix& free_transport_jacobian,
        free_vector& free_to_path_derivative) const;
};

}  // namespace detail

}  // namespace detray

#include "detray/propagator/detail/covariance_engine.ipp"