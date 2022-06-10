/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

namespace unit_constants {

// Length, native unit mm
constexpr double um = 0.001;
constexpr double mm = 1.0;
constexpr double cm = 10.0;
constexpr double m = 1000.0;

// Volume, native unit mm3
constexpr double mm3 = mm * mm * mm;
constexpr double cm3 = cm * cm * cm;

// Energy, native unit GeV
constexpr double keV = 0.000001;
constexpr double MeV = 0.001;
constexpr double GeV = 1.0;

// Atomic mass unit u
// 1u == 0.93149410242 GeV/c
constexpr double u = 0.93149410242;

// Mass
//     1eV/c² == 1.782662e-36kg
//    1GeV/c² == 1.782662e-27kg
// ->     1kg == (1/1.782662e-27)GeV/c²
// ->      1g == (1/(1e3*1.782662e-27))GeV/c²
constexpr double g = 1.0 / 1.782662e-24;
constexpr double kg = 1.0 / 1.782662e-27;

// Amount of substance, native unit mol
constexpr double mol = 1.0;
// Avogadro constant
constexpr double kAvogadro = 6.02214076e23 / unit_constants::mol;

// Magnetic field, native unit GeV/(e*mm)
constexpr double T = 0.000299792458;  // equivalent to c in appropriate SI units

}  // namespace unit_constants

}  // namespace detray
