/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/detail/density_effect_data.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <ratio>
#include <sstream>

namespace detray {

namespace detail {

/// Tag a class as 'homogeneous material', i.e. a material slab or rod.
struct homogeneous_material_tag {};

}  // namespace detail

/// Material State
enum class material_state {
    e_solid = 0,
    e_liquid = 1,
    e_gas = 2,
    e_unknown = 3
};

template <typename scalar_t, typename R = std::ratio<1, 1>>
struct material {
    using ratio = R;
    using scalar_type = scalar_t;

    constexpr material() = default;

    DETRAY_HOST_DEVICE
    constexpr material(const scalar_type x0, const scalar_type l0,
                       const scalar_type ar, const scalar_type z,
                       const scalar_type mass_rho, const material_state state)
        : m_x0(x0),
          m_l0(l0),
          m_ar(ar),
          m_z(z),
          m_mass_rho(mass_rho),
          m_state(state) {

        m_molar_rho = mass_to_molar_density(ar, mass_rho);
    }

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
          m_density(d_a, d_m, d_X0, d_X1, d_I, d_nC, d_delta0),
          m_has_density_effect_data(true) {

        m_molar_rho = mass_to_molar_density(ar, mass_rho);
    }

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr bool operator==(const material<scalar_t> &rhs) const {
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
    /// @returns the material state.
    DETRAY_HOST_DEVICE
    constexpr material_state state() const { return m_state; }
    /// @returns the molar electron density.
    DETRAY_HOST_DEVICE
    constexpr scalar_type molar_electron_density() const {
        return m_z * m_molar_rho;
    }

    /// @returns the density effect data
    DETRAY_HOST_DEVICE
    constexpr const detail::density_effect_data<scalar_type>
        &density_effect_data() const {
        return m_density;
    }

    /// @returns the (Approximated) mean excitation energy
    DETRAY_HOST_DEVICE
    scalar_type mean_excitation_energy() const {
        if (!m_has_density_effect_data) {
            // use approximative computation as defined in ATL-SOFT-PUB-2008-003
            return 16.f * unit<scalar_type>::eV *
                   math::pow(m_z, static_cast<scalar_type>(0.9));
        } else {
            return m_density.get_mean_excitation_energy();
        }
    }

    DETRAY_HOST_DEVICE
    constexpr scalar_type fraction() const {
        if constexpr (ratio::num == 0) {
            return 0.f;
        } else if constexpr (ratio::den == 0) {
            return detail::invalid_value<scalar_type>();
        } else {
            return static_cast<scalar_type>(ratio::num) /
                   static_cast<scalar_type>(ratio::den);
        }
    }

    /// @returns a string that contains the material details
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream strm;
        strm << "material: ";
        strm << " X0: " << m_x0;
        strm << " | L0: " << m_l0;
        strm << " | Z: " << m_z;

        strm << " | state: ";
        switch (m_state) {
            case material_state::e_solid: {
                strm << "solid";
                break;
            }
            case material_state::e_liquid: {
                strm << "liquid";
                break;
            }
            case material_state::e_gas: {
                strm << "gaseous";
                break;
            }
            default: {
                strm << "unknown";
                break;
            }
        };

        return strm.str();
    }

    /// @returns a string stream that prints the material details
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os, const material &mat) {
        os << mat.to_string();
        return os;
    }

    DETRAY_HOST_DEVICE
    bool has_density_effect_data() const { return m_has_density_effect_data; }

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
    scalar_type m_x0 = detail::invalid_value<scalar_type>();
    scalar_type m_l0 = detail::invalid_value<scalar_type>();
    scalar_type m_ar = 0.f;
    scalar_type m_z = 0.f;
    scalar_type m_mass_rho = 0.f;
    scalar_type m_molar_rho = 0.f;
    material_state m_state = material_state::e_unknown;
    detail::density_effect_data<scalar_type> m_density = {};
    bool m_has_density_effect_data = false;
};

// Macro for declaring the predefined materials (w/o Density effect data)
#define DETRAY_DECLARE_MATERIAL(MATERIAL_NAME, X0, L0, Ar, Z, Rho, State)   \
    template <typename scalar_t, typename R = std::ratio<1, 1>>             \
    struct MATERIAL_NAME final : public material<scalar_t, R> {             \
        using base_type = material<scalar_t, R>;                            \
        using base_type::base_type;                                         \
        DETRAY_HOST_DEVICE                                                  \
        constexpr MATERIAL_NAME() : base_type(X0, L0, Ar, Z, Rho, State) {} \
    }

// !EXPERIMENTAL!
// Macro for declaring the predefined materials (with Density effect data)
#define DETRAY_DECLARE_MATERIAL_WITH_DED(                                    \
    MATERIAL_NAME, X0, L0, Ar, Z, Rho, State, Density0, Density1, Density2,  \
    Density3, Density4, Density5, Density6)                                  \
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
