/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

namespace unit_constants {

// Length, native unit um
constexpr double um = 0.001;

// Length, native unit mm
constexpr double mm = 1.0;

// Length, convert to cm
constexpr double cm = 10.0;

// Length, convert to m
constexpr double m = 1000.0;

// Energy/mass/momentum, native unit GeV
constexpr double GeV = 1.0;

// Magnetic field, native unit GeV/(e*mm)
constexpr double T = 0.000299792458;  // equivalent to c in appropriate SI units

}  // namespace unit_constants

}  // namespace detray