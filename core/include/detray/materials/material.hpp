/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/definitions/units.hpp"

// System include(s)
#include <climits>
#include <ratio>

namespace detray {

template <typename scalar_t, typename R>
struct material {
    using ratio = R;
    using scalar_type = scalar_t;

    material(const scalar_type x0, const scalar_type l0, const scalar_type ar,
             const scalar_type z, const scalar_type mass_rho)
        : _x0(x0), _l0(l0), _ar(ar), _z(z) {

        double atomic_mass = static_cast<double>(ar) * unit_constants::u;
        _molar_rho = (mass_rho == 0)
                         ? 0
                         : static_cast<scalar_type>(
                               static_cast<double>(mass_rho) /
                               (atomic_mass * unit_constants::kAvogadro));
    }

    /// Return the radition length. Infinity in case of vacuum.
    constexpr scalar_type X0() const { return _x0; }
    /// Return the nuclear interaction length. Infinity in case of vacuum.
    constexpr scalar_type L0() const { return _l0; }
    /// Return the relative atomic mass.
    constexpr scalar_type Ar() const { return _ar; }
    /// Return the nuclear charge number.
    constexpr scalar_type Z() const { return _z; }
    /// Return the molar density.
    constexpr scalar_type molar_density() const { return _molar_rho; }
    /// Return the molar electron density.
    constexpr scalar_type molar_electron_density() const {
        return _z * _molar_rho;
    }

    protected:
    scalar_type _x0;
    scalar_type _l0;
    scalar_type _ar;
    scalar_type _z;
    scalar_type _molar_rho;
};

}  // namespace detray