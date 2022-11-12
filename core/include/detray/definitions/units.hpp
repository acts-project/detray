/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/// Unit conversion factors
template <typename scalar_t>
struct unit {

    /// Length, native unit mm
    /// @{
    static constexpr scalar_t um{0.001};
    static constexpr scalar_t mm{1.0};
    static constexpr scalar_t cm{10.0};
    static constexpr scalar_t m{1000.0};
    /// @}

    /// Volume, native unit mm3
    /// @{
    static constexpr scalar_t mm3{mm * mm * mm};
    static constexpr scalar_t cm2{cm * cm};
    static constexpr scalar_t cm3{cm * cm * cm};
    /// @}

    /// Time, native unit mm{[speed-of-light * time]{mm/s * s}}
    /// @{
    static constexpr scalar_t s{299792458000.0};
    static constexpr scalar_t fs{1e-15 * s};
    static constexpr scalar_t ps{1e-12 * s};
    static constexpr scalar_t ns{1e-9 * s};
    static constexpr scalar_t us{1e-6 * s};
    static constexpr scalar_t ms{1e-3 * s};
    static constexpr scalar_t min{60.0 * s};
    static constexpr scalar_t h{3600.0 * s};
    /// @}

    /// Energy, native unit GeV
    /// @{
    static constexpr scalar_t eV{1e-9};
    static constexpr scalar_t keV{1e-6};
    static constexpr scalar_t MeV{1e-3};
    static constexpr scalar_t GeV{1.0};
    static constexpr scalar_t TeV{1e3};
    /// @}

    /// Atomic mass unit u
    /// 1u == 0.93149410242 GeV/c
    static constexpr scalar_t u{0.93149410242};

    /// Mass
    ///     1eV/c² == 1.782662e-36kg
    ///    1GeV/c² == 1.782662e-27kg
    /// ->     1kg == (1/1.782662e-27)GeV/c²
    /// ->      1g == (1/(1e3*1.782662e-27))GeV/c²
    /// @{
    static constexpr scalar_t g{1.0 / 1.782662e-24};
    static constexpr scalar_t kg{1.0 / 1.782662e-27};
    /// @}

    /// Amount of substance, native unit mol
    static constexpr scalar_t mol{1.0};

    /// Charge, native unit e (elementary charge)
    static constexpr scalar_t e{1.0};

    /// Magnetic field, native unit GeV/(e*mm)
    static constexpr scalar_t T{
        0.000299792458};  // equivalent to c in appropriate SI units
};

/// Physical constants
template <typename scalar_t>
struct constant {

    /// Avogadro constant
    static constexpr scalar_t avogadro{6.02214076e23 / unit<scalar_t>::mol};

    /// Reduced Planck constant h/2*pi.
    ///
    /// Computed from CODATA 2018 constants to double precision.
    static constexpr scalar_t hbar{6.582119569509066e-25 * unit<scalar_t>::GeV *
                                   unit<scalar_t>::s};
};

}  // namespace detray
