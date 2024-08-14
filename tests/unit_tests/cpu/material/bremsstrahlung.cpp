/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/pdg_particle.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/test/common/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;
using algebra_t = test::algebra;

// Test class for energy loss with Bremsstrahlung
// Input tuple: < material, particle type, kinetic energy, expected output >
class EnergyLossBremsValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, pdg_particle<scalar>, scalar, scalar>> {
};

// This tests the material functionalities
TEST_P(EnergyLossBremsValidation, bremsstrahlung) {

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

    // Bremsstrahlung stopping power in MeV * cm^2 / g
    const scalar dEdx{
        I.compute_bremsstrahlung(mat, ptc, {ptc, qop}) / mat.mass_density() /
        (unit<scalar>::MeV * unit<scalar>::cm2 / unit<scalar>::g)};

    const scalar expected_dEdx = std::get<3>(GetParam());

    // We have not implemented the bremsstrahlung for the heavier charged
    // particles, which is negligible
    if (ptc.pdg_num() == electron<scalar>().pdg_num()) {
        // Check if difference is within 11% error
        EXPECT_NEAR((expected_dEdx - dEdx) / dEdx, 0.f, 0.11f);
    } else {
        EXPECT_FLOAT_EQ(static_cast<float>(dEdx), 0.f);
    }
}

// REFERENCE
//
// Atomic Data and Nuclear Data Tables Volume 4, March 1972, Pages 1-27,
//"Energy loss, range, and bremsstrahlung yield for 10-keV to 100-MeV electrons
// in various elements and chemical compounds"

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bremsstrahlung_100MeV_He, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 0.95886f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bremsstrahlung_100MeV_Al, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 3.8172f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bremsstrahlung_100MeV_Cu, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 7.2365f)));

// We have not implemented the bremsstrahlung for muons
INSTANTIATE_TEST_SUITE_P(
    detray_material_Bremsstrahlung_100MeV_Cu_muon, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), muon<scalar>(),
                                      100.0f * unit<scalar>::MeV, 0.f)));
