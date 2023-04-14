/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/concepts/algebra.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/unit_vectors.hpp"

// System include(s).
#include <random>

namespace detray {

template <CONSTRAINT(concepts::algebra) algebra_t>
struct scattering_helper {
    public:
    using matrix_operator = typename algebra_t::matrix_actor;
    using vector3 = typename algebra_t::vector3;
    using scalar_type = typename algebra_t::scalar_type;

    /// @brief Operator to scatter the direction with scattering angle
    ///
    /// @param dir  input direction
    /// @param angle  scattering angle
    /// @param generator random generator
    /// @returns the new direction from random scattering
    template <typename generator_t>
    DETRAY_HOST inline vector3 operator()(const vector3& dir,
                                          const scalar_type angle,
                                          generator_t& generator) const {

        // Generate theta and phi for random scattering
        const scalar_type r_theta{
            std::normal_distribution<scalar_type>(0.f, angle)(generator)};
        const scalar_type r_phi{std::uniform_real_distribution<scalar_type>(
            -constant<scalar_type>::pi, constant<scalar_type>::pi)(generator)};

        // xaxis of curvilinear plane
        const vector3 u = unit_vectors<vector3>().make_curvilinear_unit_u(dir);

        vector3 new_dir = axis_rotation<algebra_t>(u, r_theta)(dir);
        return axis_rotation<algebra_t>(dir, r_phi)(new_dir);
    }
};

}  // namespace detray
