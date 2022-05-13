/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <climits>

namespace detray {

/**
 * Elements Declaration
 * @note: Values from
 * https://pdg.lbl.gov/2020/AtomicNuclearProperties/index.html (Last revised 04
 * June 2020)
 */

// vacuum
DETRAY_DECLARE_MATERIAL(vacuum, std::numeric_limits<scalar>::infinity(),
                        std::numeric_limits<scalar>::infinity(), 0, 0, 0);

// H (1): Hydrogen Gas
DETRAY_DECLARE_MATERIAL(hydrogen_gas, 7.526E6 * unit_constants::mm,
                        6.209E6 * unit_constants::mm, 1.008, 1,
                        8.376E-5 * unit_constants::g /
                            (1 * unit_constants::cm3));

// H (1): Hydrogen liguid
DETRAY_DECLARE_MATERIAL(hydrogen_liquid, 8904 * unit_constants::mm,
                        7346 * unit_constants::mm, 1.008, 1,
                        0.07080 * unit_constants::g /
                            (1 * unit_constants::cm3));

// He (2): Helium Gas
DETRAY_DECLARE_MATERIAL(helium_gas, 5.671E6 * unit_constants::mm,
                        4.269E6 * unit_constants::mm, 4.003, 2,
                        1.663E-4 * unit_constants::g /
                            (1 * unit_constants::cm3));

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium, 352.8 * unit_constants::mm,
                        421.0 * unit_constants::mm, 9.012, 4.0,
                        1.848 * unit_constants::g / (1 * unit_constants::cm3));

// C (6): Carbon_gas
DETRAY_DECLARE_MATERIAL(carbon_gas, 213.5 * unit_constants::mm,
                        429.0 * unit_constants::mm, 12.01, 6,
                        2.0 * unit_constants::g / (1 * unit_constants::cm3));

// N (7): Nitrogen_gas
DETRAY_DECLARE_MATERIAL(nitrogen_gas, 3.260E+05 * unit_constants::mm,
                        7.696E+05 * unit_constants::mm, 14.007, 7,
                        1.165E-03 * unit_constants::g /
                            (1 * unit_constants::cm3));

// O (8): Oxygen Gas
DETRAY_DECLARE_MATERIAL(oxygen_gas, 2.571E+05 * unit_constants::mm,
                        6.772E+05 * unit_constants::mm, 15.999, 8,
                        1.332E-3 * unit_constants::g /
                            (1 * unit_constants::cm3));

// O (8): Oxygen liquid
DETRAY_DECLARE_MATERIAL(oxygen_liquid, 300.1 * unit_constants::mm,
                        790.3 * unit_constants::mm, 15.999, 8,
                        1.141 * unit_constants::g / (1 * unit_constants::cm3));

// Al (13)
DETRAY_DECLARE_MATERIAL(aluminium, 88.97 * unit_constants::mm,
                        397.0 * unit_constants::mm, 26.98, 13,
                        2.699 * unit_constants::g / (1 * unit_constants::cm3));

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon, 93.7 * unit_constants::mm,
                        465.2 * unit_constants::mm, 28.0855, 14.,
                        2.329 * unit_constants::g / (1 * unit_constants::cm3));

// Ar (18): Argon gas
DETRAY_DECLARE_MATERIAL(argon_gas, 1.176E+05 * unit_constants::mm,
                        7.204E+05 * unit_constants::mm, 39.948, 18.,
                        1.662E-03 * unit_constants::g /
                            (1 * unit_constants::cm3));

// W (74)
DETRAY_DECLARE_MATERIAL(tungsten, 3.504 * unit_constants::mm,
                        99.46 * unit_constants::mm, 183.84, 74,
                        19.3 * unit_constants::g / (1 * unit_constants::cm3));

// Au (79)
DETRAY_DECLARE_MATERIAL(gold, 3.344 * unit_constants::mm,
                        101.6 * unit_constants::mm, 196.97, 79,
                        19.32 * unit_constants::g / (1 * unit_constants::cm3));

/**
 * Elements Declaration for ACTS Generic detector
 * @note: Values taken from BuildGenericDetector.hpp in ACTS
 */

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium_tml, 352.8 * unit_constants::mm,
                        407.0 * unit_constants::mm, 9.012, 4.0,
                        1.848 * unit_constants::g / (1 * unit_constants::cm3));

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon_tml, 95.7 * unit_constants::mm,
                        465.2 * unit_constants::mm, 28.03, 14.,
                        2.32 * unit_constants::g / (1 * unit_constants::cm3));

/**
 * Mixtures or Compounds
 */

// Air (dry, 1 atm)
// @note:
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/air_dry_1_atm.html
// @note: Ar from Wikipedia
DETRAY_DECLARE_MATERIAL(air, 3.039E+05 * unit_constants::mm,
                        7.477E+05 * unit_constants::mm, 28.97, 14.46,
                        1.205E-03 * unit_constants::g /
                            (1 * unit_constants::cm3));

// (CH3)2CHCH3 Gas
// @note: (X0, L0, mass_rho) from https://pdg.lbl.gov/2005/reviews/atomicrpp.pdf
// @note: Ar from Wikipedia
// @note: Z was caculated by simply summing the number of atoms. Surprisingly
// it seems the right value because Z/A is 0.58496, which is the same with <Z/A>
// in the pdg refernce
DETRAY_DECLARE_MATERIAL(isobutane, 169300 * unit_constants::mm,
                        288.3 * unit_constants::mm, 58.124, 34.,
                        2.67 * unit_constants::g / (1 * unit_constants::cm3));

// C3H8 Gas
// @note: (X0, L0, mass_rho) from
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/propane.html
// @note: Ar from Wikipedia
DETRAY_DECLARE_MATERIAL(propane, 2.429E+05 * unit_constants::mm,
                        4.106E+05 * unit_constants::mm, 44.097, 26,
                        1.868E-03 * unit_constants::g /
                            (1 * unit_constants::cm3));

}  // namespace detray