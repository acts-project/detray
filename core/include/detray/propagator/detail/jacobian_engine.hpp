/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracks/detail/track_helper.hpp"

namespace detray {

namespace detail {

// Jacobian engine assuming "Planar Surface with Cartesian coordinate"
// NOTE: We may inherit jacobian engine to address different
// types of surface, if required
template <typename transform3_t>
struct jacobian_engine {

    /// Transformation matching this struct
    using transform3_type = transform3_t;
    /// scalar_type
    using scalar_type = typename transform3_type::scalar_type;
    /// Vector in 3D space
    using vector3 = typename transform3_type::vector3;
    /// Matrix operator
    using matrix_operator = typename transform3_type::matrix_actor;
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;
    using size_type = typename transform3_type::size_type;
    /// 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    /// Shorthand vector types related to track parameters.
    using bound_vector = matrix_type<e_bound_size, 1>;
    using free_vector = matrix_type<e_free_size, 1>;
    /// Mapping matrix
    using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
    using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
    using free_to_path_matrix = matrix_type<1, e_free_size>;

    /** Function to get the jacobian for bound to free coordinate transform
     *
     * @param trf3 is the transform matrix
     * @param bound_vec is the input bound vector
     * @returns bound to free jacobian
     */
    DETRAY_HOST_DEVICE inline bound_to_free_matrix bound_to_free_coordinate(
        const transform3_type& trf3, const bound_vector& bound_vec) const;

    /** Function to get the jacobian for free to bound coordinate transform
     *
     * @param trf3 is the transform matrix
     * @param free_vec is the input free vector
     * @returns free to bound jacobian
     */
    DETRAY_HOST_DEVICE inline free_to_bound_matrix free_to_bound_coordinate(
        const transform3_type& trf3, const free_vector& free_vec) const;

    /** Function to calculate the path correction term for transport jacobian,
     * which is caused by the geometry constraint.
     *
     * @param trf3 is the transform matrix
     * @param free_vec is the input free vector
     * @returns free to path matrix
     */
    DETRAY_HOST_DEVICE inline free_to_path_matrix free_to_path_correction(
        const transform3_type& trf3, const free_vector& free_vec) const;
};

}  // namespace detail

}  // namespace detray

#include "detray/propagator/detail/jacobian_engine.ipp"