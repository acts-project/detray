/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/detail/density_effect_data.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <limits>

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
                        std::numeric_limits<scalar>::infinity(), 0.f, 0.f, 0.f,
                        material_state::e_gas, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
                        0.f);

// H₂ (1): Hydrogen Gas
DETRAY_DECLARE_MATERIAL(hydrogen_gas, 7.526E3f * unit<scalar>::m,
                        6.209E3f * unit<scalar>::m, 2.016f, 2.f,
                        static_cast<scalar>(8.376E-5 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_gas, 0.1409f, 5.7273f, 1.8639f,
                        3.2718f, 19.2f, 9.5834f, 0.f);

// H₂ (1): Hydrogen Liquid
DETRAY_DECLARE_MATERIAL(hydrogen_liquid, 8.904f * unit<scalar>::m,
                        7.346f * unit<scalar>::m, 2.016f, 2.f,
                        static_cast<scalar>(0.07080f * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_liquid, 0.1348f, 5.6249f, 0.4400f,
                        1.8856f, 21.8f, 3.0977f, 0.f);

// He (2): Helium Gas
DETRAY_DECLARE_MATERIAL(helium_gas, 5.671E3f * unit<scalar>::m,
                        4.269E3f * unit<scalar>::m, 4.003f, 2.f,
                        static_cast<scalar>(1.663E-4 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_gas, 0.1344f, 5.8347f, 2.2017f,
                        3.6122f, 41.8f, 11.1393f, 0.f);

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium, 352.8f * unit<scalar>::mm,
                        421.0f * unit<scalar>::mm, 9.012f, 4.f,
                        static_cast<scalar>(1.848 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_solid, 0.8039f, 2.4339f, 0.0592f,
                        1.6922f, 63.7f, 2.7847f, 0.14f);

// C (6): Carbon (amorphous)
DETRAY_DECLARE_MATERIAL(
    carbon_gas, 213.5f * unit<scalar>::mm, 429.0f * unit<scalar>::mm, 12.01f,
    6.f, static_cast<scalar>(2.0 * unit<double>::g / unit<double>::cm3),
    material_state::e_gas, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

// N₂ (7): Nitrogen Gas
DETRAY_DECLARE_MATERIAL(nitrogen_gas, 3.260E+02f * unit<scalar>::m,
                        7.696E+02f * unit<scalar>::m, 28.014f, 14.f,
                        static_cast<scalar>(1.165E-03 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_gas, 0.1535f, 3.2125f, 1.7378f,
                        4.1323f, 82.0f, 10.5400f, 0.f);

// O₂ (8): Oxygen Gas
DETRAY_DECLARE_MATERIAL(oxygen_gas, 2.571E+02f * unit<scalar>::m,
                        6.772E+02f * unit<scalar>::m, 31.998f, 16.f,
                        static_cast<scalar>(1.332E-3 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_gas, 0.1178f, 3.2913f, 1.7541f,
                        4.3213f, 95.0f, 10.7004f, 0.f);

// O₂ (8): Oxygen liquid
DETRAY_DECLARE_MATERIAL(oxygen_liquid, 300.1f * unit<scalar>::mm,
                        790.3f * unit<scalar>::mm, 31.998f, 16.f,
                        static_cast<scalar>(1.141 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_liquid, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
                        0.f);

// Al (13)
DETRAY_DECLARE_MATERIAL(aluminium, 88.97f * unit<scalar>::mm,
                        397.0f * unit<scalar>::mm, 26.98f, 13.f,
                        static_cast<scalar>(2.699 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_solid, 0.0802f, 3.6345f, 0.1708f,
                        3.0127f, 166.f, 4.2395f, 0.12f);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon, 93.7f * unit<scalar>::mm,
                        465.2f * unit<scalar>::mm, 28.0855f, 14.f,
                        static_cast<scalar>(2.329 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_solid, 0.1492f, 3.2546f, 0.2015f,
                        2.8716f, 173.0f, 4.4355f, 0.14f);

// Ar (18): Argon gas
DETRAY_DECLARE_MATERIAL(argon_gas, 1.176E+02f * unit<scalar>::m,
                        7.204E+02f * unit<scalar>::m, 39.948f, 18.f,
                        static_cast<scalar>(1.662E-03 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_gas, 0.1971f, 2.9618f, 1.7635f,
                        4.4855f, 188.f, 11.9480f, 0.f);

// W (74)
DETRAY_DECLARE_MATERIAL(tungsten, 3.504f * unit<scalar>::mm,
                        99.46f * unit<scalar>::mm, 183.84f, 74.f,
                        static_cast<scalar>(19.3 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_solid, 0.1551f, 2.8447f, 0.2167f,
                        3.4960f, 727.0f, 5.4059f, 0.14f);

// Au (79)
DETRAY_DECLARE_MATERIAL(gold, 3.344f * unit<scalar>::mm,
                        101.6f * unit<scalar>::mm, 196.97f, 79.f,
                        static_cast<scalar>(19.32 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_solid, 0.0976f, 3.1101f, 0.2021f,
                        3.6979f, 790.f, 5.5747f, 0.14f);

/**
 * Elements Declaration for ACTS Generic detector
 * @note: Values taken from BuildGenericDetector.hpp in ACTS
 */

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium_tml, 352.8f * unit<scalar>::mm,
                        407.f * unit<scalar>::mm, 9.012f, 4.f,
                        static_cast<scalar>(1.848 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_solid, 0.8039f, 2.4339f, 0.0592f,
                        1.6922f, 63.7f, 2.7847f, 0.14f);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon_tml, 95.7f * unit<scalar>::mm,
                        465.2f * unit<scalar>::mm, 28.03f, 14.f,
                        static_cast<scalar>(2.32 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_solid, 0.1492f, 3.2546f, 0.2015f,
                        2.8716f, 173.f, 4.4355f, 0.14f);

/**
 * Mixtures or Compounds
 */

// Air (dry, 1 atm)
// @note:
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/air_dry_1_atm.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Molar_mass)
DETRAY_DECLARE_MATERIAL(air, 3.039E+02f * unit<scalar>::m,
                        7.477E+02f * unit<scalar>::m, 28.97f, 14.46f,
                        static_cast<scalar>(1.205E-03 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_gas, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
                        0.f);

// (CH3)2CHCH3 Gas
// @note: (X0, L0, mass_rho) from https://pdg.lbl.gov/2005/reviews/atomicrpp.pdf
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Isobutane)
// @note: Z was caculated by simply summing the number of atoms. Surprisingly
// it seems the right value because Z/A is 0.58496, which is the same with <Z/A>
// in the pdg refernce
DETRAY_DECLARE_MATERIAL(
    isobutane, 1693E+02f * unit<scalar>::mm, 288.3f * unit<scalar>::mm, 58.124f,
    34.f, static_cast<scalar>(2.67 * unit<double>::g / unit<double>::cm3),
    material_state::e_gas, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

// C3H8 Gas
// @note: (X0, L0, mass_rho) from
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/propane.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Propane)
DETRAY_DECLARE_MATERIAL(propane, 2.429E+02f * unit<scalar>::m,
                        4.106E+02f * unit<scalar>::m, 44.097f, 26.f,
                        static_cast<scalar>(1.868E-03 * unit<double>::g /
                                            unit<double>::cm3),
                        material_state::e_gas, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
                        0.f);

}  // namespace detray