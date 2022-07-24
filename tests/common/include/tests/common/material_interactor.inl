/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

/// Detray include(s)
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/interactor.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/predefined_materials.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using point2 = __plugin::point2<scalar>;
using point3 = __plugin::point3<scalar>;

// Test class for MUON energy loss with Bethe function
// Input tuple: < material / energy / expected output from
// https://pdg.lbl.gov/2022/AtomicNuclearProperties for Muon dEdX and range >
class EnergyLossBetheValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, scalar, scalar>> {};

// This tests the material functionalities
TEST_P(EnergyLossBetheValidation, bethe_energy_loss) {

    // Interaction object
    interactor<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // H2 liquid with a unit thickness
    material_slab<scalar> slab(std::get<0>(GetParam()), 1 * unit_constants::cm);

    // muon
    const int pdg = pdg_particle::eMuon;

    // mass
    const scalar m = 105.7 * unit_constants::MeV;

    // qOverP
    const scalar qOverP = -1. / std::get<1>(GetParam());

    // Bethe Stopping power in MeV * cm^2 / g
    const scalar dEdx =
        I.compute_energy_loss_bethe(is, slab, pdg, m, qOverP) /
        slab.path_segment(is) / slab.get_material().mass_density() /
        (unit_constants::MeV * unit_constants::cm2 / unit_constants::g);

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - dEdx) / dEdx < 0.05);
}

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      0.1003 * unit_constants::GeV, 6.539)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      1.101 * unit_constants::GeV, 4.182)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      10.11 * unit_constants::GeV, 4.777)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      100.1 * unit_constants::GeV, 5.305)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      0.1003 * unit_constants::GeV, 3.082)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      1.101 * unit_constants::GeV, 2.133)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      10.11 * unit_constants::GeV, 2.768)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      100.1 * unit_constants::GeV, 3.188)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      0.1003 * unit_constants::GeV, 2.533)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      1.101 * unit_constants::GeV, 1.744)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      10.11 * unit_constants::GeV, 2.097)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      100.1 * unit_constants::GeV, 2.360)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      0.1003 * unit_constants::GeV, 2.608)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      1.101 * unit_constants::GeV, 1.803)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      10.11 * unit_constants::GeV, 2.177)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      100.1 * unit_constants::GeV, 2.451)));

// Test class for MUON energy loss with Landau function
// Input tuple: < material / energy / expected energy loss  / expected fwhm  >
class EnergyLossLandauValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, scalar, scalar, scalar>> {};

TEST_P(EnergyLossLandauValidation, landau_energy_loss) {

    // Interaction object
    interactor<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // H2 liquid with a unit thickness
    material_slab<scalar> slab(std::get<0>(GetParam()),
                               0.17 * unit_constants::cm);

    // muon
    const int pdg = pdg_particle::eMuon;

    // mass
    const scalar m = 105.7 * unit_constants::MeV;

    // qOverP
    const scalar qOverP = -1. / std::get<1>(GetParam());

    // Landau Energy loss in MeV
    const scalar dE = I.compute_energy_loss_landau(is, slab, pdg, m, qOverP) /
                      (unit_constants::MeV);

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - dE) / dE < 0.05);

    // Landau Energy loss Fluctuation
    const scalar fwhm =
        I.compute_energy_loss_landau_fwhm(is, slab, pdg, m, qOverP) /
        (unit_constants::MeV);

    // Check if difference is within 10% error
    EXPECT_TRUE(std::abs(std::get<3>(GetParam()) - fwhm) / fwhm < 0.1);
}

// Expected output from Fig 33.7 in RPP2018
INSTANTIATE_TEST_SUITE_P(
    Landau_10GeV_Silicon, EnergyLossLandauValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      10. * unit_constants::GeV, 0.525, 0.13)));