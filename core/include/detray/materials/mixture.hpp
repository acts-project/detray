/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/materials/material.hpp"
#include "detray/utils/ratio.hpp"

// System include(s)
#include <tuple>

namespace detray {

// Compile-time material mixture type. The summation of ratios should be eqaul
// to one
template <typename scalar_t, typename... material_types>
struct mixture
    : public material<scalar_t, typename ratio_sum<
                                    typename material_types::ratio...>::ratio> {
    public:
    using ratio = typename ratio_sum<typename material_types::ratio...>::ratio;

    using scalar_type = scalar_t;
    using base_type = material<scalar_type, ratio>;
    using base_type::base_type;

    static_assert(is_ratio_one_v<ratio>,
                  "Sumation of ratios should be equal to 1");

    // Constructor
    mixture() {
        // Compute effective relative atomic mass
        // Zeff = Ar0 * ratio0 + Ar1 * ratio1 + ...
        auto sum_Ar = [](material_types... M) -> decltype(auto) {
            return ((M.Ar() * M.fraction()) + ...);
        };

        this->m_ar = std::apply(sum_Ar, std::tuple<material_types...>());

        // Compute effective atomic number
        // Zeff = Z0 * ratio0 + Z1 * ratio1 + ...
        auto sum_Z = [](material_types... M) -> decltype(auto) {
            return ((M.Z() * M.fraction()) + ...);
        };

        this->m_z = std::apply(sum_Z, std::tuple<material_types...>());

        // Get averaged mass density
        auto sum_rho = [](material_types... M) -> decltype(auto) {
            return ((M.mass_density() * M.fraction()) + ...);
        };

        this->m_mass_rho = std::apply(sum_rho, std::tuple<material_types...>());

        // Compute effective radiation length (X0)
        // reference:
        // https://cds.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
        // X_avg = W_avg / Sum_i[ W_i / X_i ],
        // where:
        // W_i is mass density of i_th component
        // W_avg is the averaged mass density
        auto sum_rho_over_X0 = [](material_types... M) -> decltype(auto) {
            return ((M.fraction() / M.X0()) + ...);
        };
        this->m_x0 =
            1. / std::apply(sum_rho_over_X0, std::tuple<material_types...>());

        // Compute effective nuclear radiation length
        // Follow the same equation of effective X0
        auto sum_rho_over_L0 = [](material_types... M) -> decltype(auto) {
            return ((M.fraction() / M.L0()) + ...);
        };

        this->m_l0 =
            1. / std::apply(sum_rho_over_L0, std::tuple<material_types...>());

        // Compute molar density
        this->m_molar_rho =
            this->mass_to_molar_density(this->m_ar, this->m_mass_rho);
    }
};

}  // namespace detray