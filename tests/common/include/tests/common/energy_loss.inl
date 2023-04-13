/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

using namespace detray;

// Project include(s).
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/predefined_materials.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// Test class for MUON energy loss with Bethe function
// Input tuple: < material / energy / expected output from
// https://pdg.lbl.gov/2022/AtomicNuclearProperties for Muon dEdX and range >
class EnergyLossBetheValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, scalar, scalar>> {};

// This tests the material functionalities
TEST_P(EnergyLossBetheValidation, bethe_energy_loss) {

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // H2 liquid with a unit thickness
    material_slab<scalar> slab(std::get<0>(GetParam()), 1.f * unit<scalar>::cm);

    // muon
    constexpr int pdg = pdg_particle::eMuon;

    // mass
    constexpr scalar m{105.7f * unit<scalar>::MeV};

    // qOverP
    const scalar qOverP{-1.f / std::get<1>(GetParam())};

    // Bethe Stopping power in MeV * cm^2 / g
    const scalar dEdx{
        I.compute_energy_loss_bethe(is, slab, pdg, m, qOverP, -1.f) /
        slab.path_segment(is) / slab.get_material().mass_density() /
        (unit<scalar>::MeV * unit<scalar>::cm2 / unit<scalar>::g)};

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - dEdx) / dEdx < 0.05f);
}

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 6.539f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      1.101f * unit<scalar>::GeV, 4.182f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      10.11f * unit<scalar>::GeV, 4.777f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      100.1f * unit<scalar>::GeV, 5.305f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 3.082f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      1.101f * unit<scalar>::GeV, 2.133f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      10.11f * unit<scalar>::GeV, 2.768f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      100.1f * unit<scalar>::GeV, 3.188f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.533f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.744f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      10.11f * unit<scalar>::GeV, 2.097f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.360f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.608f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.803f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 2.177f)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.451f)));

// Declare Silicon without density effect data
constexpr const material<scalar> Si_approx(
    93.7f * unit<scalar>::mm, 465.2f * unit<scalar>::mm, 28.0855f, 14.f,
    static_cast<scalar>(2.329 * unit<double>::g / unit<double>::cm3),
    material_state::e_solid, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

INSTANTIATE_TEST_SUITE_P(Bethe_0p1GeV_Si_approx, EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             Si_approx, 0.1003f * unit<scalar>::GeV, 2.608f)));

INSTANTIATE_TEST_SUITE_P(Bethe_1GeV_Si_approx, EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             Si_approx, 1.101f * unit<scalar>::GeV, 1.803f)));

INSTANTIATE_TEST_SUITE_P(Bethe_10GeV_Si_approx, EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             Si_approx, 10.11f * unit<scalar>::GeV, 2.177f)));

INSTANTIATE_TEST_SUITE_P(Bethe_100GeV_Si_approx, EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             Si_approx, 100.1f * unit<scalar>::GeV, 2.451f)));

// Test class for MUON energy loss with Landau function
// Input tuple: < material / energy / expected energy loss  / expected fwhm  >
class EnergyLossLandauValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, scalar, scalar, scalar>> {};

TEST_P(EnergyLossLandauValidation, landau_energy_loss) {

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // H2 liquid with a unit thickness
    material_slab<scalar> slab(std::get<0>(GetParam()),
                               0.17f * unit<scalar>::cm);

    // muon
    constexpr int pdg{pdg_particle::eMuon};

    // mass
    constexpr scalar m{105.7f * unit<scalar>::MeV};

    // qOverP
    const scalar qOverP{-1.f / std::get<1>(GetParam())};

    // Landau Energy loss in MeV
    const scalar dE{
        I.compute_energy_loss_landau(is, slab, pdg, m, qOverP, -1.f) /
        unit<scalar>::MeV};

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - dE) / dE < 0.05f);

    // Landau Energy loss Fluctuation
    const scalar fwhm{
        I.compute_energy_loss_landau_fwhm(is, slab, pdg, m, qOverP, -1.f) /
        unit<scalar>::MeV};

    // Check if difference is within 10% error
    EXPECT_TRUE(std::abs(std::get<3>(GetParam()) - fwhm) / fwhm < 0.1f);
}

// Expected output from Fig 33.7 in RPP2018
INSTANTIATE_TEST_SUITE_P(Landau_10GeV_Silicon, EnergyLossLandauValidation,
                         ::testing::Values(std::make_tuple(
                             silicon<scalar>(), 10.f * unit<scalar>::GeV,
                             0.525f, 0.13f)));