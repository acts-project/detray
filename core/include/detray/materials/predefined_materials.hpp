/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/materials/detail/density_effect_data.hpp"
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

// Density effect data taken from
// https://pdg.lbl.gov/2022/AtomicNuclearProperties for Muon dEdX and range

// Vacuum
DETRAY_DECLARE_MATERIAL(vacuum, std::numeric_limits<scalar>::infinity(),
                        std::numeric_limits<scalar>::infinity(), 0, 0, 0,
                        material_state::e_gas, 0, 0, 0, 0, 0, 0, 0);

// H₂ (1): Hydrogen Gas
DETRAY_DECLARE_MATERIAL(hydrogen_gas, 7.526E6 * unit<scalar>::mm,
                        6.209E6 * unit<scalar>::mm, 1.008 * 2, 1 * 2,
                        8.376E-5 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0.1409, 5.7273, 1.8639, 3.2718,
                        19.2, 9.5834, 0.00);

// H₂ (1): Hydrogen Liquid
DETRAY_DECLARE_MATERIAL(hydrogen_liquid, 8904 * unit<scalar>::mm,
                        7346 * unit<scalar>::mm, 1.008 * 2, 1 * 2,
                        0.07080 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_liquid, 0.1348, 5.6249, 0.4400,
                        1.8856, 21.8, 3.0977, 0.00);

// He (2): Helium Gas
DETRAY_DECLARE_MATERIAL(helium_gas, 5.671E6 * unit<scalar>::mm,
                        4.269E6 * unit<scalar>::mm, 4.003, 2,
                        1.663E-4 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0.1344, 5.8347, 2.2017, 3.6122,
                        41.8, 11.1393, 0.00);

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium, 352.8 * unit<scalar>::mm,
                        421.0 * unit<scalar>::mm, 9.012, 4.0,
                        1.848 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_solid, 0.8039, 2.4339, 0.0592, 1.6922,
                        63.7, 2.7847, 0.14);

// C (6): Carbon (amorphous)
DETRAY_DECLARE_MATERIAL(carbon_gas, 213.5 * unit<scalar>::mm,
                        429.0 * unit<scalar>::mm, 12.01, 6,
                        2.0 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0, 0, 0, 0, 0, 0, 0);

// N₂ (7): Nitrogen Gas
DETRAY_DECLARE_MATERIAL(nitrogen_gas, 3.260E+05 * unit<scalar>::mm,
                        7.696E+05 * unit<scalar>::mm, 14.007 * 2, 7 * 2,
                        1.165E-03 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0.1535, 3.2125, 1.7378, 4.1323,
                        82.0, 10.5400, 0.00);

// O₂ (8): Oxygen Gas
DETRAY_DECLARE_MATERIAL(oxygen_gas, 2.571E+05 * unit<scalar>::mm,
                        6.772E+05 * unit<scalar>::mm, 15.999 * 2, 8 * 2,
                        1.332E-3 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0.1178, 3.2913, 1.7541, 4.3213,
                        95.0, 10.7004, 0.00);

// O₂ (8): Oxygen liquid
DETRAY_DECLARE_MATERIAL(oxygen_liquid, 300.1 * unit<scalar>::mm,
                        790.3 * unit<scalar>::mm, 15.999 * 2, 8 * 2,
                        1.141 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_liquid, 0, 0, 0, 0, 0, 0, 0);

// Al (13)
DETRAY_DECLARE_MATERIAL(aluminium, 88.97 * unit<scalar>::mm,
                        397.0 * unit<scalar>::mm, 26.98, 13,
                        2.699 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_solid, 0.0802, 3.6345, 0.1708, 3.0127,
                        166.0, 4.2395, 0.12);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon, 93.7 * unit<scalar>::mm,
                        465.2 * unit<scalar>::mm, 28.0855, 14.,
                        2.329 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_solid, 0.1492, 3.2546, 0.2015, 2.8716,
                        173.0, 4.4355, 0.14);

// Ar (18): Argon gas
DETRAY_DECLARE_MATERIAL(argon_gas, 1.176E+05 * unit<scalar>::mm,
                        7.204E+05 * unit<scalar>::mm, 39.948, 18.,
                        1.662E-03 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0.1971, 2.9618, 1.7635, 4.4855,
                        188.0, 11.9480, 0.00);

// W (74)
DETRAY_DECLARE_MATERIAL(tungsten, 3.504 * unit<scalar>::mm,
                        99.46 * unit<scalar>::mm, 183.84, 74,
                        19.3 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_solid, 0.1551, 2.8447, 0.2167, 3.4960,
                        727.0, 5.4059, 0.14);

// Au (79)
DETRAY_DECLARE_MATERIAL(gold, 3.344 * unit<scalar>::mm,
                        101.6 * unit<scalar>::mm, 196.97, 79,
                        19.32 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_solid, 0.0976, 3.1101, 0.2021, 3.6979,
                        790.0, 5.5747, 0.14);

/**
 * Elements Declaration for ACTS Generic detector
 * @note: Values taken from BuildGenericDetector.hpp in ACTS
 */

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium_tml, 352.8 * unit<scalar>::mm,
                        407.0 * unit<scalar>::mm, 9.012, 4.0,
                        1.848 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_solid, 0.8039, 2.4339, 0.0592, 1.6922,
                        63.7, 2.7847, 0.14);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon_tml, 95.7 * unit<scalar>::mm,
                        465.2 * unit<scalar>::mm, 28.03, 14.,
                        2.32 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_solid, 0.1492, 3.2546, 0.2015, 2.8716,
                        173.0, 4.4355, 0.14);

/**
 * Mixtures or Compounds
 */

// Air (dry, 1 atm)
// @note:
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/air_dry_1_atm.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Molar_mass)
DETRAY_DECLARE_MATERIAL(air, 3.039E+05 * unit<scalar>::mm,
                        7.477E+05 * unit<scalar>::mm, 28.97, 14.46,
                        1.205E-03 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0, 0, 0, 0, 0, 0, 0);

// (CH3)2CHCH3 Gas
// @note: (X0, L0, mass_rho) from https://pdg.lbl.gov/2005/reviews/atomicrpp.pdf
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Isobutane)
// @note: Z was caculated by simply summing the number of atoms. Surprisingly
// it seems the right value because Z/A is 0.58496, which is the same with <Z/A>
// in the pdg refernce
DETRAY_DECLARE_MATERIAL(isobutane, 169300 * unit<scalar>::mm,
                        288.3 * unit<scalar>::mm, 58.124, 34.,
                        2.67 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0, 0, 0, 0, 0, 0, 0);

// C3H8 Gas
// @note: (X0, L0, mass_rho) from
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/propane.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Propane)
DETRAY_DECLARE_MATERIAL(propane, 2.429E+05 * unit<scalar>::mm,
                        4.106E+05 * unit<scalar>::mm, 44.097, 26,
                        1.868E-03 * unit<scalar>::g / (1 * unit<scalar>::cm3),
                        material_state::e_gas, 0, 0, 0, 0, 0, 0, 0);

}  // namespace detray