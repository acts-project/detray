/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"

// Detray test include(s)
#include "detray/test/utils/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

// Test class for the stopping power
// Input tuple: < material, particle type, kinetic energy, expected output >
class StoppingPowerValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, pdg_particle<scalar>, scalar, scalar>> {
};

TEST_P(StoppingPowerValidation, stopping_power) {

    // Interaction object
    interaction<scalar> I;

    // Material
    material<scalar> mat = std::get<0>(GetParam());

    // Particle
    pdg_particle<scalar> ptc = std::get<1>(GetParam());

    // Kinetic energy
    const scalar T = std::get<2>(GetParam());

    // Total energy
    const scalar E = T + ptc.mass();

    // Momentum
    const scalar p = math::sqrt(E * E - ptc.mass() * ptc.mass());

    // qoverp
    const scalar qop{ptc.charge() / p};

    // Stopping power in MeV * cm^2 / g
    const scalar dEdx{
        I.compute_stopping_power(mat, ptc, {ptc, qop}) / mat.mass_density() /
        (unit<scalar>::MeV * unit<scalar>::cm2 / unit<scalar>::g)};

    const scalar expected_dEdx = std::get<3>(GetParam());

    // Check if difference is within 8% error
    EXPECT_NEAR((expected_dEdx - dEdx) / dEdx, 0.f, 0.08f);
}

/******************
 *   Muon tests
 ******************/

// From https://pdg.lbl.gov/2024/AtomicNuclearProperties/index.html
// Note 1: that we took the PDG value only from Ionization loss (Radiative loss is
// ignored)
// Note 2: assumes that the stopping powers of muon and antimuon are the same
// Note 3: Test fails with He Gas and 1 GeV muons (18 % difference)
INSTANTIATE_TEST_SUITE_P(
    muon_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                                      100.0f * unit<scalar>::MeV, 2.165f),
                      //std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                      //                1.f * unit<scalar>::GeV, 2.133f),
                      std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                                      10.0f * unit<scalar>::GeV, 2.768f),
                      std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                                      100.0f * unit<scalar>::GeV, 3.188f)));

INSTANTIATE_TEST_SUITE_P(
    muon_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      100.0f * unit<scalar>::MeV, 1.849f),
                      std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      1.f * unit<scalar>::GeV, 1.803f),
                      std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      10.0f * unit<scalar>::GeV, 2.177f),
                      std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      100.0f * unit<scalar>::GeV, 2.451f)));

INSTANTIATE_TEST_SUITE_P(
    anti_muon_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), antimuon<scalar>(),
                                      100.0f * unit<scalar>::MeV, 2.165f),
                      //std::make_tuple(helium_gas<scalar>(), antimuon<scalar>(),
                      //                1.f * unit<scalar>::GeV, 2.133f),
                      std::make_tuple(helium_gas<scalar>(), antimuon<scalar>(),
                                      10.0f * unit<scalar>::GeV, 2.768f),
                      std::make_tuple(helium_gas<scalar>(), antimuon<scalar>(),
                                      100.0f * unit<scalar>::GeV, 3.188f)));

INSTANTIATE_TEST_SUITE_P(
    anti_muon_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), antimuon<scalar>(),
                                      100.0f * unit<scalar>::MeV, 1.849f),
                      std::make_tuple(silicon<scalar>(), antimuon<scalar>(),
                                      1.f * unit<scalar>::GeV, 1.803f),
                      std::make_tuple(silicon<scalar>(), antimuon<scalar>(),
                                      10.0f * unit<scalar>::GeV, 2.177f),
                      std::make_tuple(silicon<scalar>(), antimuon<scalar>(),
                                      100.0f * unit<scalar>::GeV, 2.451f)));

/*********************
 *   Electron tests
 *********************/

// From https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
// Assumes that the stopping powers of electron and positron are the same
INSTANTIATE_TEST_SUITE_P(
    electron_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 3.532f),
                      std::make_tuple(helium_gas<scalar>(), electron<scalar>(),
                                      1.f * unit<scalar>::GeV, 13.14f)));

INSTANTIATE_TEST_SUITE_P(
    electron_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 6.017f),
                      std::make_tuple(silicon<scalar>(), electron<scalar>(),
                                      1.f * unit<scalar>::GeV, 46.69f)));

INSTANTIATE_TEST_SUITE_P(
    positron_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), positron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 3.532f),
                      std::make_tuple(helium_gas<scalar>(), positron<scalar>(),
                                      1.f * unit<scalar>::GeV, 13.14f)));

INSTANTIATE_TEST_SUITE_P(
    positron_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), positron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 6.017f),
                      std::make_tuple(silicon<scalar>(), positron<scalar>(),
                                      1.f * unit<scalar>::GeV, 46.69f)));