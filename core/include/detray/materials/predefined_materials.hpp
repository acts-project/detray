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

// Vacuum
DETRAY_DECLARE_MATERIAL(vacuum, std::numeric_limits<scalar>::infinity(),
                        std::numeric_limits<scalar>::infinity(), 0, 0, 0,
                        material_state::e_gas, detail::null_density);

// H₂ (1): Hydrogen Gas
DETRAY_DECLARE_MATERIAL(hydrogen_gas, 7.526E6 * unit_constants::mm,
                        6.209E6 * unit_constants::mm, 1.008 * 2, 1 * 2,
                        8.376E-5 * unit_constants::g /
                            (1 * unit_constants::cm3),
                        material_state::e_gas, detail::H2_gas_density);

// H₂ (1): Hydrogen Liquid
DETRAY_DECLARE_MATERIAL(hydrogen_liquid, 8904 * unit_constants::mm,
                        7346 * unit_constants::mm, 1.008 * 2, 1 * 2,
                        0.07080 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_liquid, detail::H2_liquid_density);

// He (2): Helium Gas
DETRAY_DECLARE_MATERIAL(helium_gas, 5.671E6 * unit_constants::mm,
                        4.269E6 * unit_constants::mm, 4.003, 2,
                        1.663E-4 * unit_constants::g /
                            (1 * unit_constants::cm3),
                        material_state::e_gas, detail::He_gas_density);

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium, 352.8 * unit_constants::mm,
                        421.0 * unit_constants::mm, 9.012, 4.0,
                        1.848 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_solid, detail::Be_density);

// C (6): Carbon (amorphous)
DETRAY_DECLARE_MATERIAL(carbon_gas, 213.5 * unit_constants::mm,
                        429.0 * unit_constants::mm, 12.01, 6,
                        2.0 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_gas, detail::null_density);

// N₂ (7): Nitrogen Gas
DETRAY_DECLARE_MATERIAL(nitrogen_gas, 3.260E+05 * unit_constants::mm,
                        7.696E+05 * unit_constants::mm, 14.007 * 2, 7 * 2,
                        1.165E-03 * unit_constants::g /
                            (1 * unit_constants::cm3),
                        material_state::e_gas, detail::N2_gas_density);

// O₂ (8): Oxygen Gas
DETRAY_DECLARE_MATERIAL(oxygen_gas, 2.571E+05 * unit_constants::mm,
                        6.772E+05 * unit_constants::mm, 15.999 * 2, 8 * 2,
                        1.332E-3 * unit_constants::g /
                            (1 * unit_constants::cm3),
                        material_state::e_gas, detail::O2_gas_density);

// O₂ (8): Oxygen liquid
DETRAY_DECLARE_MATERIAL(oxygen_liquid, 300.1 * unit_constants::mm,
                        790.3 * unit_constants::mm, 15.999 * 2, 8 * 2,
                        1.141 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_liquid, detail::null_density);

// Al (13)
DETRAY_DECLARE_MATERIAL(aluminium, 88.97 * unit_constants::mm,
                        397.0 * unit_constants::mm, 26.98, 13,
                        2.699 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_solid, detail::Al_density);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon, 93.7 * unit_constants::mm,
                        465.2 * unit_constants::mm, 28.0855, 14.,
                        2.329 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_solid, detail::Si_density);

// Ar (18): Argon gas
DETRAY_DECLARE_MATERIAL(argon_gas, 1.176E+05 * unit_constants::mm,
                        7.204E+05 * unit_constants::mm, 39.948, 18.,
                        1.662E-03 * unit_constants::g /
                            (1 * unit_constants::cm3),
                        material_state::e_gas, detail::Ar_gas_density);

// W (74)
DETRAY_DECLARE_MATERIAL(tungsten, 3.504 * unit_constants::mm,
                        99.46 * unit_constants::mm, 183.84, 74,
                        19.3 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_solid, detail::W_density);

// Au (79)
DETRAY_DECLARE_MATERIAL(gold, 3.344 * unit_constants::mm,
                        101.6 * unit_constants::mm, 196.97, 79,
                        19.32 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_solid, detail::Au_density);

/**
 * Elements Declaration for ACTS Generic detector
 * @note: Values taken from BuildGenericDetector.hpp in ACTS
 */

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium_tml, 352.8 * unit_constants::mm,
                        407.0 * unit_constants::mm, 9.012, 4.0,
                        1.848 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_solid, detail::Be_density);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon_tml, 95.7 * unit_constants::mm,
                        465.2 * unit_constants::mm, 28.03, 14.,
                        2.32 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_solid, detail::Si_density);

/**
 * Mixtures or Compounds
 */

// Air (dry, 1 atm)
// @note:
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/air_dry_1_atm.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Molar_mass)
DETRAY_DECLARE_MATERIAL(air, 3.039E+05 * unit_constants::mm,
                        7.477E+05 * unit_constants::mm, 28.97, 14.46,
                        1.205E-03 * unit_constants::g /
                            (1 * unit_constants::cm3),
                        material_state::e_gas, detail::null_density);

// (CH3)2CHCH3 Gas
// @note: (X0, L0, mass_rho) from https://pdg.lbl.gov/2005/reviews/atomicrpp.pdf
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Isobutane)
// @note: Z was caculated by simply summing the number of atoms. Surprisingly
// it seems the right value because Z/A is 0.58496, which is the same with <Z/A>
// in the pdg refernce
DETRAY_DECLARE_MATERIAL(isobutane, 169300 * unit_constants::mm,
                        288.3 * unit_constants::mm, 58.124, 34.,
                        2.67 * unit_constants::g / (1 * unit_constants::cm3),
                        material_state::e_gas, detail::null_density);

// C3H8 Gas
// @note: (X0, L0, mass_rho) from
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/propane.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Propane)
DETRAY_DECLARE_MATERIAL(propane, 2.429E+05 * unit_constants::mm,
                        4.106E+05 * unit_constants::mm, 44.097, 26,
                        1.868E-03 * unit_constants::g /
                            (1 * unit_constants::cm3),
                        material_state::e_gas, detail::null_density);

}  // namespace detray