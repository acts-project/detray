/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/detail/density_effect_data.hpp"

// System include(s)
#include <limits>
#include <ratio>

namespace detray {

/// Material State
enum class material_state { e_solid = 0, e_liquid = 1, e_gas = 2 };

template <typename scalar_t, typename R = std::ratio<1, 1>>
struct material {
    using ratio = R;
    using scalar_type = scalar_t;

    constexpr material() = default;

    DETRAY_HOST_DEVICE
    constexpr material(const scalar_type x0, const scalar_type l0,
                       const scalar_type ar, const scalar_type z,
                       const scalar_type mass_rho, const material_state state,
                       const scalar_type d_a, const scalar_type d_m,
                       const scalar_type d_X0, const scalar_type d_X1,
                       const scalar_type d_I, const scalar_type d_nC,
                       const scalar_type d_delta0)
        : m_x0(x0),
          m_l0(l0),
          m_ar(ar),
          m_z(z),
          m_mass_rho(mass_rho),
          m_state(state),
          m_density(d_a, d_m, d_X0, d_X1, d_I, d_nC, d_delta0) {

        m_molar_rho = mass_to_molar_density(ar, mass_rho);
    }

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr bool operator==(const material<scalar_t>& rhs) const {
        return (m_x0 == rhs.X0() && m_l0 == rhs.L0() && m_ar == rhs.Ar() &&
                m_z == rhs.Z());
    }

    /// @returns the radition length. Infinity in case of vacuum.
    DETRAY_HOST_DEVICE
    constexpr scalar_type X0() const { return m_x0; }
    /// @returns the nuclear interaction length. Infinity in case of vacuum.
    DETRAY_HOST_DEVICE
    constexpr scalar_type L0() const { return m_l0; }
    /// @returns the relative atomic mass.
    DETRAY_HOST_DEVICE
    constexpr scalar_type Ar() const { return m_ar; }
    /// @returns the nuclear charge number.
    DETRAY_HOST_DEVICE
    constexpr scalar_type Z() const { return m_z; }
    /// @returns the mass density.
    DETRAY_HOST_DEVICE
    constexpr scalar_type mass_density() const { return m_mass_rho; }
    /// @returns the molar density.
    DETRAY_HOST_DEVICE
    constexpr scalar_type molar_density() const { return m_molar_rho; }
    /// @returns the molar electron density.
    DETRAY_HOST_DEVICE
    constexpr scalar_type molar_electron_density() const {
        return m_z * m_molar_rho;
    }
    /// @returns the density effect data
    DETRAY_HOST_DEVICE
    constexpr detail::density_effect_data<scalar_type> density_effect_data()
        const {
        return m_density;
    }

    /// @returns the (Approximated) mean excitation energy
    DETRAY_HOST_DEVICE
    scalar_type mean_excitation_energy() const {
        // use approximative computation as defined in ATL-SOFT-PUB-2008-003
        if (m_density == detail::density_effect_data<scalar_type>{}) {
            return 16.f * unit<scalar_type>::eV *
                   math_ns::pow(m_z, static_cast<scalar_type>(0.9));
        } else {
            return m_density.get_mean_excitation_energy();
        }
    }

    DETRAY_HOST_DEVICE
    constexpr scalar_type fraction() const {
        if constexpr (ratio::num == 0) {
            return 0.f;
        } else if constexpr (ratio::den == 0) {
            return std::numeric_limits<scalar_type>::infinity();
        } else {
            return static_cast<scalar_type>(ratio::num) /
                   static_cast<scalar_type>(ratio::den);
        }
    }

    protected:
    DETRAY_HOST_DEVICE
    constexpr scalar_type mass_to_molar_density(double ar, double mass_rho) {
        if (mass_rho == 0.) {
            return 0.f;
        }

        const double molar_mass{ar * unit<double>::u *
                                constant<double>::avogadro};

        return static_cast<scalar_type>(mass_rho / molar_mass);
    }

    // Material properties
    scalar_type m_x0 = std::numeric_limits<scalar_type>::infinity();
    scalar_type m_l0 = std::numeric_limits<scalar_type>::infinity();
    scalar_type m_ar = 0.f;
    scalar_type m_z = 0.f;
    scalar_type m_mass_rho = 0.f;
    scalar_type m_molar_rho = 0.f;
    material_state m_state = material_state::e_solid;
    detail::density_effect_data<scalar_type> m_density = {};
};

// Macro for declaring the predefined materials (with Density effect data)
#define DETRAY_DECLARE_MATERIAL(MATERIAL_NAME, X0, L0, Ar, Z, Rho, State,    \
                                Density0, Density1, Density2, Density3,      \
                                Density4, Density5, Density6)                \
    template <typename scalar_t, typename R = std::ratio<1, 1>>              \
    struct MATERIAL_NAME final : public material<scalar_t, R> {              \
        using base_type = material<scalar_t, R>;                             \
        using base_type::base_type;                                          \
        DETRAY_HOST_DEVICE                                                   \
        constexpr MATERIAL_NAME()                                            \
            : base_type(X0, L0, Ar, Z, Rho, State, Density0, Density1,       \
                        Density2, Density3, Density4, Density5, Density6) {} \
    }

}  // namespace detray