/** Detray library, part of the ACTS project
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/propagator/detail/jacobian.hpp"

namespace detray::detail {

/// @brief Specialization for 2D concentric cylindrical frames
template <concepts::algebra algebra_t>
struct jacobian<concentric_cylindrical2D<algebra_t>> {

    /// @name Type definitions for the struct
    /// @{
    using coordinate_frame = concentric_cylindrical2D<algebra_t>;

    using algebra_type = algebra_t;
    using transform3_type = dtransform3D<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;

    using bound_to_free_matrix_type = bound_to_free_matrix<algebra_t>;
    using free_to_bound_matrix_type = free_to_bound_matrix<algebra_t>;
    using free_to_path_matrix_type = free_to_path_matrix<algebra_t>;

    DETRAY_HOST_DEVICE static constexpr free_to_path_matrix_type
    path_derivative(const transform3_type & /*trf*/, const point3_type &pos,
                    const vector3_type &dir, const vector3_type & /*dtds*/) {

        using mask_t = mask<concentric_cylinder2D, algebra_type, std::uint8_t>;
        constexpr mask_t dummy_mask{};

        const transform3_type identity{};

        const point2_type local =
            coordinate_frame::global_to_local(identity, pos, {});
        const vector3_type normal =
            coordinate_frame::normal(identity, local, dummy_mask);

        const vector3_type pos_term =
            (-1.f / vector::dot(normal, dir)) * normal;

        auto derivative{matrix::zero<free_to_path_matrix_type>()};
        getter::element(derivative, 0u, e_free_pos0) = pos_term[0];
        getter::element(derivative, 0u, e_free_pos1) = pos_term[1];
        getter::element(derivative, 0u, e_free_pos2) = pos_term[2];

        return derivative;
    }

    DETRAY_HOST_DEVICE
    static constexpr void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix_type &bound_to_free_jacobian,
        const transform3_type & /*trf*/, const point3_type &pos,
        const vector3_type & /*dir*/) {

        const scalar_type r{vector::perp(pos)};
        const scalar_type phi{vector::phi(pos)};

        getter::element(bound_to_free_jacobian, e_free_pos0, 0u) =
            -r * math::sin(phi);
        getter::element(bound_to_free_jacobian, e_free_pos1, 0u) =
            r * math::cos(phi);
        getter::element(bound_to_free_jacobian, e_free_pos2, e_bound_loc1) =
            1.f;
    }

    DETRAY_HOST_DEVICE
    static constexpr void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix_type &free_to_bound_jacobian,
        const transform3_type & /*trf*/, const point3_type &pos,
        const vector3_type & /*dir*/) {

        const scalar_type r_inv{1.f / vector::perp(pos)};
        const scalar_type phi{vector::phi(pos)};

        getter::element(free_to_bound_jacobian, 0u, e_free_pos0) =
            -r_inv * math::sin(phi);
        getter::element(free_to_bound_jacobian, 0u, e_free_pos1) =
            r_inv * math::cos(phi);
        getter::element(free_to_bound_jacobian, e_bound_loc1, e_free_pos2) =
            1.f;
    }

    DETRAY_HOST_DEVICE
    static constexpr void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix_type & /*bound_to_free_jacobian*/,
        const transform3_type & /*trf3*/, const point3_type & /*pos*/,
        const vector3_type & /*dir*/) {
        // Do nothing
    }
};

}  // namespace detray::detail
