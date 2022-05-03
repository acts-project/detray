/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <climits>

namespace detray {

template <typename scalar_t>
struct material {
    using scalar_type = scalar_t;

    material(const scalar_type x0, const scalar_type l0, const scalar_type ar,
             const scalar_type z, const scalar_type molar_rho)
        : _x0(x0), _l0(l0), _ar(ar), _z(z), _molar_rho(molar_rho) {}

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

    private:
    scalar_type _x0 = std::numeric_limits<scalar_type>::infinity();
    scalar_type _l0 = std::numeric_limits<scalar_type>::infinity();
    scalar_type _ar = 0.0f;
    scalar_type _z = 0.0f;
    scalar_type _molar_rho = 0.0f;
};

}  // namespace detray