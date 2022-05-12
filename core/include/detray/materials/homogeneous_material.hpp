/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/materials/material.hpp"

namespace detray {

template <typename scalar_t>
struct homogeneous_material : public material<scalar_t, std::ratio<1, 1>> {
    using base_type = material<scalar_t, std::ratio<1, 1>>;
    using base_type::base_type;
    using scalar_type = scalar_t;

    homogeneous_material(const scalar_type x0, const scalar_type l0,
                         const scalar_type ar, const scalar_type z,
                         const scalar_type mass_rho)
        : base_type(x0, l0, ar, z, mass_rho) {}
};

#define DETRAY_DECLARE_MATERIAL(MATERIAL_NAME, X0, L0, Ar, Z, Rho)       \
    template <typename scalar_t>                                         \
    struct MATERIAL_NAME final : public homogeneous_material<scalar_t> { \
        using base_type = homogeneous_material<scalar_t>;                \
        using base_type::base_type;                                      \
        MATERIAL_NAME() : base_type(X0, L0, Ar, Z, Rho) {}               \
    };

}  // namespace detray