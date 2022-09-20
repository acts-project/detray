/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/detail/density_effect_data.hpp"

// System include(s)
#include <climits>
#include <ratio>

namespace detray {

// Material State
enum class material_state { e_solid = 0, e_liquid = 1, e_gas = 2 };

template <typename scalar_t, typename R = std::ratio<1, 1>>
struct material {
    using ratio = R;
    using scalar_type = scalar_t;

    material() = default;

    DETRAY_HOST_DEVICE
    material(const scalar_type x0, const scalar_type l0, const scalar_type ar,
             const scalar_type z, const scalar_type mass_rho,
             const material_state state, const scalar_type d_a,
             const scalar_type d_m, const scalar_type d_X0,
             const scalar_type d_X1, const scalar_type d_I,
             const scalar_type d_nC, const scalar_type d_delta0)
        : m_x0(x0),
          m_l0(l0),
          m_ar(ar),
          m_z(z),
          m_mass_rho(mass_rho),
          m_state(state),
          m_density(d_a, d_m, d_X0, d_X1, d_I, d_nC, d_delta0) {

        m_molar_rho = mass_to_molar_density(ar, mass_rho);
    }

    /** Equality operator
     *
     * @param rhs is the right hand side to be compared to
     */
    DETRAY_HOST_DEVICE bool operator==(const material<scalar_t> &rhs) const {
        return (m_x0 == rhs.X0() && m_l0 == rhs.L0() && m_ar == rhs.Ar() &&
                m_z == rhs.Z());
    }

    /// Return the radition length. Infinity in case of vacuum.
    DETRAY_HOST_DEVICE
    constexpr scalar_type X0() const { return m_x0; }
    /// Return the nuclear interaction length. Infinity in case of vacuum.
    DETRAY_HOST_DEVICE
    constexpr scalar_type L0() const { return m_l0; }
    /// Return the relative atomic mass.
    DETRAY_HOST_DEVICE
    constexpr scalar_type Ar() const { return m_ar; }
    /// Return the nuclear charge number.
    DETRAY_HOST_DEVICE
    constexpr scalar_type Z() const { return m_z; }
    /// Return the mass density.
    DETRAY_HOST_DEVICE
    constexpr scalar_type mass_density() const { return m_mass_rho; }
    /// Return the molar density.
    DETRAY_HOST_DEVICE
    constexpr scalar_type molar_density() const { return m_molar_rho; }
    /// Return the molar electron density.
    DETRAY_HOST_DEVICE
    constexpr scalar_type molar_electron_density() const {
        return m_z * m_molar_rho;
    }
    /// Return the density effect data
    DETRAY_HOST_DEVICE
    constexpr detail::density_effect_data<scalar_type> density_effect_data()
        const {
        return m_density;
    }

    /// Return the (Approximated) mean excitation energy
    DETRAY_HOST_DEVICE
    scalar_type mean_excitation_energy() const {
        // use approximative computation as defined in ATL-SOFT-PUB-2008-003
        if (m_density ==
            detail::density_effect_data<scalar_type>(0, 0, 0, 0, 0, 0, 0)) {
            return scalar_type(16 * unit_constants::eV) *
                   std::pow(m_z, scalar_type(0.9));
        } else {
            return m_density.get_mean_excitation_energy();
        }
    }

    DETRAY_HOST_DEVICE
    constexpr scalar_type fraction() const {
        return ratio::num / static_cast<scalar_type>(ratio::den);
    }

    protected:
    DETRAY_HOST_DEVICE
    scalar_type mass_to_molar_density(scalar_type ar, scalar_type mass_rho) {
        if (mass_rho == 0) {
            return 0;
        }

        double atomic_mass = static_cast<double>(ar) * unit_constants::u;

        return static_cast<scalar_type>(
            static_cast<double>(mass_rho) /
            (atomic_mass * unit_constants::kAvogadro));
    }

    // Material properties
    scalar_type m_x0 = std::numeric_limits<scalar_type>::infinity();
    scalar_type m_l0 = std::numeric_limits<scalar_type>::infinity();
    scalar_type m_ar = 0;
    scalar_type m_z = 0;
    scalar_type m_mass_rho = 0;
    scalar_type m_molar_rho = 0;
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
        MATERIAL_NAME()                                                      \
            : base_type(X0, L0, Ar, Z, Rho, State, Density0, Density1,       \
                        Density2, Density3, Density4, Density5, Density6) {} \
    }

}  // namespace detray