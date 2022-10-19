/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/// Symbolic values for commonly used PDG particle numbers.
enum pdg_particle : int32_t {
    eInvalid = 0,
    eElectron = 11,
    eAntiElectron = -eElectron,
    ePositron = -eElectron,
    eMuon = 13,
    eAntiMuon = -eMuon,
    eTau = 15,
    eAntiTau = -eTau,
    eGamma = 22,
    ePionZero = 111,
    ePionPlus = 211,
    ePionMinus = -ePionPlus,
    eNeutron = 2112,
    eAntiNeutron = -eNeutron,
    eProton = 2212,
    eAntiProton = -eProton,
};

}  // namespace detray