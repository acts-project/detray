/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/test/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

GTEST_TEST(derivation_test, beta) {

    // mass
    constexpr scalar m{105.7f * unit<scalar>::MeV};

    // charge
    const scalar q = -1.f;

    // Displacement for numerical differentiaion
    const scalar h = 0.001f;

    // Iterate from 1 GeV to 10 GeV
    for (unsigned int i = 1u; i < 10; i++) {
        const scalar p = static_cast<scalar>(i) * detray::unit<scalar>::GeV;
        const scalar qop = q / p;

        detail::relativistic_quantities rq(m, qop, q);
        detail::relativistic_quantities rq1(m, qop + h, q);
        detail::relativistic_quantities rq2(m, qop - h, q);

        const scalar numerical = (rq1.m_beta - rq2.m_beta) / (2.f * h);

        const scalar evaluated = rq.derive_beta();

        EXPECT_NEAR(numerical, evaluated, numerical * 0.01);
    }
}

// Test class for the derivation of bethe stopping power
class DerivationStoppingPowerValidation
    : public ::testing::TestWithParam<std::tuple<material<scalar>>> {};

TEST_P(DerivationStoppingPowerValidation, derivation_of_stopping_power) {

    // Interaction object
    interaction<scalar> I;

    // Material
    const auto mat = std::get<0>(GetParam());

    // muon
    constexpr int pdg = pdg_particle::eMuon;

    // mass
    constexpr scalar m{105.7f * unit<scalar>::MeV};

    // charge
    const scalar q = -1.f;

    // Displacement for numerical differentiaion
    const scalar h = 0.0001f;

    // Iterate from 2 GeV to 100 GeV
    for (unsigned int i = 2u; i < 100; i++) {
        const scalar p = static_cast<scalar>(i) * detray::unit<scalar>::GeV;
        const scalar qop = q / p;

        const scalar g1 = I.compute_stopping_power(mat, pdg, m, qop + h, q);
        const scalar g2 = I.compute_stopping_power(mat, pdg, m, qop - h, q);

        const scalar numerical = (g1 - g2) / (2.f * h);

        const scalar evaluated = I.derive_stopping_power(mat, pdg, m, qop, q);

        EXPECT_NEAR(numerical, evaluated, numerical * 0.01);
    }
}

INSTANTIATE_TEST_SUITE_P(StoppingPowerDerivationSilicon,
                         DerivationStoppingPowerValidation,
                         ::testing::Values(silicon<scalar>()));
