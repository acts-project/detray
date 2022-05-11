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
#include "detray/materials/material_composition.hpp"

// System include(s)
#include <climits>

namespace detray {

/**
 * @note: Values taken from pdg.lbl.gov/2020/AtomicNuclearProperties/index.html
 * (Last revised 04 June 2020 )
 */

// vacuum
DETRAY_DECLARE_MATERIAL(vacuum, std::numeric_limits<scalar>::infinity(),
                        std::numeric_limits<scalar>::infinity(), 0, 0, 0);

// H (1): Hydrogen Gas
DETRAY_DECLARE_MATERIAL(hydrogen, 7.526E6 * unit_constants::mm,
                        6.209E6 * unit_constants::mm, 1.008, 1,
                        8.376E-5 * unit_constants::g /
                            (1 * unit_constants::cm3));

// He (2): Helium Gas
DETRAY_DECLARE_MATERIAL(helium, 5.671E6 * unit_constants::mm,
                        4.269E6 * unit_constants::mm, 4.003, 2,
                        1.663E-4 * unit_constants::g /
                            (1 * unit_constants::cm3));

// Be (4)
DETRAY_DECLARE_MATERIAL(beryllium, 352.8 * unit_constants::mm,
                        421.0 * unit_constants::mm, 9.012, 4.0,
                        1.848 * unit_constants::g / (1 * unit_constants::cm3));

// C (6)
DETRAY_DECLARE_MATERIAL(carbon, 213.5 * unit_constants::mm,
                        429.0 * unit_constants::mm, 12.01, 6,
                        2.0 * unit_constants::g / (1 * unit_constants::cm3));

// Al (13)
DETRAY_DECLARE_MATERIAL(aluminium, 88.97 * unit_constants::mm,
                        397.0 * unit_constants::mm, 26.98, 13,
                        2.699 * unit_constants::g / (1 * unit_constants::cm3));

// Si (14)
DETRAY_DECLARE_MATERIAL(silicon, 93.7 * unit_constants::mm,
                        465.2 * unit_constants::mm, 28.0855, 14.,
                        2.329 * unit_constants::g / (1 * unit_constants::cm3));

// W (74)
DETRAY_DECLARE_MATERIAL(tungsten, 3.504 * unit_constants::mm,
                        99.46 * unit_constants::mm, 183.84, 74,
                        19.3 * unit_constants::g / (1 * unit_constants::cm3));

// Au (79)
DETRAY_DECLARE_MATERIAL(gold, 3.344 * unit_constants::mm,
                        101.6 * unit_constants::mm, 196.97, 79,
                        19.32 * unit_constants::g / (1 * unit_constants::cm3));

/**
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

}  // namespace detray