/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/simulation/landau_distribution.hpp"
#include "detray/simulation/random_scatterer.hpp"
#include "detray/test/types.hpp"
#include "detray/utils/statistics.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <random>

using namespace detray;
using transform3 = test::transform3;

using sf_desc_t = surface_descriptor<>;

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
    intersection2D<sf_desc_t> is;
    is.cos_incidence_angle = 1.f;

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

    const scalar expected_dEdx = std::get<2>(GetParam());

    // Check if difference is within 5% error
    EXPECT_NEAR((expected_dEdx - dEdx) / dEdx, 0.f, 0.05f);
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

/*
//@ NOTE: Test fails with He Gas and 10 GeV muons (18 % difference)
INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      1.101f * unit<scalar>::GeV, 2.133f)));
*/

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

// Test class for MUON energy loss with Landau function
// Input tuple: < material / energy / expected energy loss  / expected fwhm  >
class EnergyLossLandauValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, scalar, scalar, scalar>> {};

TEST_P(EnergyLossLandauValidation, landau_energy_loss) {

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    intersection2D<sf_desc_t> is;
    is.cos_incidence_angle = 1.f;

    // Material
    const auto mat = std::get<0>(GetParam());

    // Thickness
    const scalar thickness = 0.17f * unit<scalar>::cm;

    // Material slab with a unit thickness
    material_slab<scalar> slab(mat, thickness);

    // muon
    constexpr int pdg{pdg_particle::eMuon};

    // mass
    constexpr scalar m{105.7f * unit<scalar>::MeV};

    // Energy
    const scalar E = std::get<1>(GetParam());

    // p
    const scalar p = std::sqrt(E * E - m * m);

    // q
    const scalar q = -1.f;

    // qOverP
    const scalar qOverP{q / p};

    // Landau Energy loss in MeV
    const scalar Landau_MeV{
        I.compute_energy_loss_landau(is, slab, pdg, m, qOverP, q) /
        unit<scalar>::MeV};

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - Landau_MeV) / Landau_MeV <
                0.05f);

    // Landau Energy loss Fluctuation
    const scalar fwhm_MeV{
        I.compute_energy_loss_landau_fwhm(is, slab, pdg, m, qOverP, q) /
        unit<scalar>::MeV};

    // Check if difference is within 10% error
    EXPECT_TRUE(std::abs(std::get<3>(GetParam()) - fwhm_MeV) / fwhm_MeV < 0.1f);
}

// Expected output from Fig 33.7 in RPP2018
INSTANTIATE_TEST_SUITE_P(Landau_10GeV_Silicon, EnergyLossLandauValidation,
                         ::testing::Values(std::make_tuple(
                             silicon<scalar>(), 10.f * unit<scalar>::GeV,
                             0.525f, 0.13f)));

// Input tuple: < energy >
class LandauDistributionValidation
    : public ::testing::TestWithParam<std::tuple<scalar>> {};

TEST_P(LandauDistributionValidation, landau_distribution) {

    // Random generator
    std::random_device rd{};
    std::mt19937_64 generator{rd()};
    generator.seed(0u);

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    intersection2D<sf_desc_t> is;
    is.cos_incidence_angle = 1.f;

    // Material
    const auto mat = silicon<scalar>();

    // Thickness
    const auto thickness = 1.f * unit<scalar>::mm;

    // Material slab with a unit thickness
    material_slab<scalar> slab(mat, thickness);

    // muon
    constexpr int pdg{pdg_particle::eMuon};

    // mass
    constexpr scalar m{105.7f * unit<scalar>::MeV};

    // Energy
    const scalar E = std::get<0>(GetParam());

    // p
    const scalar p = std::sqrt(E * E - m * m);

    // q
    const scalar q = -1.f;

    // qOverP
    const scalar qOverP{q / p};

    // Bethe energy loss
    const scalar dE{I.compute_energy_loss_bethe(is, slab, pdg, m, qOverP, q)};

    // Landau Energy loss
    const scalar mpv{I.compute_energy_loss_landau(is, slab, pdg, m, qOverP, q)};

    // Landau Energy loss Sigma
    const scalar sigma{
        I.compute_energy_loss_landau_sigma(is, slab, pdg, m, qOverP, q)};

    // Landau Sampling
    landau_distribution<scalar> landau;

    std::vector<scalar> samples;
    std::size_t n_samples{10000u};
    for (std::size_t i = 0u; i < n_samples; i++) {
        const scalar sa = landau(generator, mpv, sigma);
        samples.push_back(sa);
    }

    // Mean energy loss from Landau samples
    const scalar sample_mean = statistics::mean(samples);

    // Make sure that the difference is within 50% (...)
    // @note: We should not expect a good consistency between landau
    // distribution samples and bethe bloch function because the arguments
    // for landau_distribution are not the most probable value and sigma
    // (sigma is not even defined in Landau distribution). We might need to put
    // a smart scaling factor
    EXPECT_NEAR((dE - sample_mean) / dE, 0.f, 0.5f);

    // Test attenuate function
    std::vector<scalar> energies;

    for (std::size_t i = 0u; i < n_samples; i++) {
        const scalar new_p = random_scatterer<transform3>().attenuate(
            mpv, sigma, m, p, generator);
        ASSERT_TRUE(new_p < p);

        const scalar new_E = std::sqrt(new_p * new_p + m * m);
        energies.push_back(new_E);
    }

    const scalar E_mean = statistics::mean(energies);
    const scalar E_expected = E - dE;
    EXPECT_TRUE(E_expected < E);

    // Make sure that the difference is within 65%
    EXPECT_NEAR((E_mean - E_expected) / dE, 0.f, 0.65f);
}

INSTANTIATE_TEST_SUITE_P(
    Landau, LandauDistributionValidation,
    ::testing::Values(std::make_tuple(1.f * unit<scalar>::GeV),
                      std::make_tuple(10.f * unit<scalar>::GeV),
                      std::make_tuple(100.f * unit<scalar>::GeV)));
