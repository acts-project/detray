/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <limits>

namespace detray {

/**
 * Elements Declaration
 * @note: Values from
 * https://pdg.lbl.gov/2020/AtomicNuclearProperties/index.html (Last revised 04
 * June 2020)
 *
 * @NOTE: For diatomic molecules (e.g. H₂ and N₂), the atomic mass (A) and
 * charge number (Z) are doubled
 */
// Vacuum
DETRAY_DECLARE_MATERIAL(vacuum, std::numeric_limits<scalar_t>::max(),
                        std::numeric_limits<scalar_t>::max(), 0., 0., 0.,
                        material_state::e_unknown);

// H₂ (1): Hydrogen Gas
DETRAY_DECLARE_MATERIAL(hydrogen_gas, 7.526E3f * unit<scalar_t>::m,
                        6.209E3f * unit<scalar_t>::m, 2. * 1.008, 2. * 1.,
                        static_cast<scalar_t>(8.376E-5 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// H₂ (1): Hydrogen Liquid
DETRAY_DECLARE_MATERIAL(hydrogen_liquid, 8.904 * unit<scalar_t>::m,
                        7.346 * unit<scalar_t>::m, 2. * 1.008, 2. * 1.,
                        static_cast<scalar_t>(0.07080 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_liquid);

// He (2): Helium Gas
DETRAY_DECLARE_MATERIAL(helium_gas, 5.671E3f * unit<scalar_t>::m,
                        4.269E3f * unit<scalar_t>::m, 4.003, 2.,
                        static_cast<scalar_t>(1.663E-4 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium, 352.8 * unit<scalar_t>::mm,
                        421.0 * unit<scalar_t>::mm, 9.012, 4.,
                        static_cast<scalar_t>(1.848 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

// C (6): Carbon (amorphous)
DETRAY_DECLARE_MATERIAL(carbon_gas, 213.5 * unit<scalar_t>::mm,
                        429.0 * unit<scalar_t>::mm, 12.01, 6.,
                        static_cast<scalar_t>(2.0 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// N₂ (7): Nitrogen Gas
DETRAY_DECLARE_MATERIAL(nitrogen_gas, 3.260E+02 * unit<scalar_t>::m,
                        7.696E+02 * unit<scalar_t>::m, 2. * 14.007, 2. * 7.,
                        static_cast<scalar_t>(1.165E-03 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// O₂ (8): Oxygen Gas
DETRAY_DECLARE_MATERIAL(oxygen_gas, 2.571E+02 * unit<scalar_t>::m,
                        6.772E+02 * unit<scalar_t>::m, 2. * 15.999, 2. * 8.,
                        static_cast<scalar_t>(1.332E-3 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// O₂ (8): Oxygen liquid
DETRAY_DECLARE_MATERIAL(oxygen_liquid, 300.1 * unit<scalar_t>::mm,
                        790.3 * unit<scalar_t>::mm, 2. * 15.999, 2. * 8.,
                        static_cast<scalar_t>(1.141 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_liquid);

// Al (13)
DETRAY_DECLARE_MATERIAL(aluminium, 88.97 * unit<scalar_t>::mm,
                        397.0 * unit<scalar_t>::mm, 26.98, 13.,
                        static_cast<scalar_t>(2.699 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon, 93.7 * unit<scalar_t>::mm,
                        465.2 * unit<scalar_t>::mm, 28.0855, 14.,
                        static_cast<scalar_t>(2.329 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

// Si (14) with density effect data
DETRAY_DECLARE_MATERIAL_WITH_DED(silicon_with_ded, 93.7 * unit<scalar_t>::mm,
                                 465.2 * unit<scalar_t>::mm, 28.0855, 14.,
                                 static_cast<scalar_t>(2.329 * unit<double>::g /
                                                       unit<double>::cm3),
                                 material_state::e_solid, 0.1492, 3.2546,
                                 0.2015, 2.8716, 173.0, 4.4355, 0.14);

// Ar (18): Argon gas
DETRAY_DECLARE_MATERIAL(argon_gas, 1.176E+02 * unit<scalar_t>::m,
                        7.204E+02 * unit<scalar_t>::m, 39.948, 18.,
                        static_cast<scalar_t>(1.662E-03 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// Ar (18): Argon liquid
DETRAY_DECLARE_MATERIAL(argon_liquid, 14. * unit<scalar_t>::cm,
                        85.77 * unit<scalar_t>::cm, 39.948, 18.,
                        static_cast<scalar_t>(1.396 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_liquid);

// Fe (26)
DETRAY_DECLARE_MATERIAL(iron, 1.757 * unit<scalar_t>::cm,
                        16.77 * unit<scalar_t>::cm, 55.845, 26.,
                        static_cast<scalar_t>(7.874 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

// Fe (26) with density effect data
DETRAY_DECLARE_MATERIAL_WITH_DED(iron_with_ded, 1.757 * unit<scalar_t>::cm,
                                 16.77 * unit<scalar_t>::cm, 55.845, 26.,
                                 static_cast<scalar_t>(7.874 * unit<double>::g /
                                                       unit<double>::cm3),
                                 material_state::e_solid, 0.14680, 2.9632,
                                 -0.0012, 3.1531, 286.0, 4.2911, 0.12);

// Copper (29)
DETRAY_DECLARE_MATERIAL(copper, 1.436 * unit<scalar_t>::cm,
                        15.32 * unit<scalar_t>::cm, 63.546, 29.,
                        static_cast<scalar_t>(8.960 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

// Copper (29) with density effect data
DETRAY_DECLARE_MATERIAL_WITH_DED(copper_with_ded, 1.436 * unit<scalar_t>::cm,
                                 15.32 * unit<scalar_t>::cm, 63.546, 29.,
                                 static_cast<scalar_t>(8.960 * unit<double>::g /
                                                       unit<double>::cm3),
                                 material_state::e_solid, 0.14339, 2.9044,
                                 -0.0254, 3.2792, 322.0, 4.4190, 0.08);

// W (74)
DETRAY_DECLARE_MATERIAL(tungsten, 3.504 * unit<scalar_t>::mm,
                        99.46 * unit<scalar_t>::mm, 183.84, 74.,
                        static_cast<scalar_t>(19.3 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

// Au (79)
DETRAY_DECLARE_MATERIAL(gold, 3.344 * unit<scalar_t>::mm,
                        101.6 * unit<scalar_t>::mm, 196.97, 79.,
                        static_cast<scalar_t>(19.32 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

/**
 * Elements Declaration for ACTS Generic detector
 * @note: Values taken from BuildGenericDetector.hpp in ACTS
 */

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium_tml, 352.8 * unit<scalar_t>::mm,
                        407. * unit<scalar_t>::mm, 9.012, 4.,
                        static_cast<scalar_t>(1.848 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon_tml, 95.7 * unit<scalar_t>::mm,
                        465.2 * unit<scalar_t>::mm, 28.03, 14.,
                        static_cast<scalar_t>(2.32 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

/**
 * Mixtures or Compounds
 */

// Air (dry, 1 atm)
// @note:
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/air_dry_1_atm.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Molar_mass)
DETRAY_DECLARE_MATERIAL(air, 3.039E+02 * unit<scalar_t>::m,
                        7.477E+02 * unit<scalar_t>::m, 28.97, 14.46,
                        static_cast<scalar_t>(1.205E-03 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// (CH3)2CHCH3 Gas
// @note: (X0, L0, mass_rho) from https://pdg.lbl.gov/2005/reviews/atomicrpp.pdf
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Isobutane)
// @note: Z was caculated by simply summing the number of atoms. Surprisingly
// it seems the right value because Z/A is 0.58496, which is the same with <Z/A>
// in the pdg refernce
DETRAY_DECLARE_MATERIAL(isobutane, 1693E+02 * unit<scalar_t>::mm,
                        288.3 * unit<scalar_t>::mm, 58.124, 34.,
                        static_cast<scalar_t>(2.67 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// C3H8 Gas
// @note: (X0, L0, mass_rho) from
// https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/propane.html
// @note: Ar from Wikipedia (https://en.wikipedia.org/wiki/Propane)
DETRAY_DECLARE_MATERIAL(propane, 2.429E+02 * unit<scalar_t>::m,
                        4.106E+02 * unit<scalar_t>::m, 44.097, 26.,
                        static_cast<scalar_t>(1.868E-03 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_gas);

// Cesium Iodide (CsI)
// https://pdg.lbl.gov/2023/AtomicNuclearProperties/HTML/cesium_iodide_CsI.html
DETRAY_DECLARE_MATERIAL(cesium_iodide, 1.86 * unit<scalar_t>::cm,
                        38.04 * unit<scalar_t>::cm, 259.81, 108.,
                        static_cast<scalar_t>(4.510 * unit<double>::g /
                                              unit<double>::cm3),
                        material_state::e_solid);

DETRAY_DECLARE_MATERIAL_WITH_DED(cesium_iodide_with_ded,
                                 1.86 * unit<scalar_t>::cm,
                                 38.04 * unit<scalar_t>::cm, 259.81, 108.,
                                 static_cast<scalar_t>(4.510 * unit<double>::g /
                                                       unit<double>::cm3),
                                 material_state::e_solid, 0.25381, 2.6657,
                                 0.0395, 3.3353, 553.1, 6.2807, 0.00);

}  // namespace detray
