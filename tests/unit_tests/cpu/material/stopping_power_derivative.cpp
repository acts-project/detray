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

using namespace detray;

GTEST_TEST(detray_material, derivative_test_beta2) {

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

        const scalar numerical = (rq1.m_beta2 - rq2.m_beta2) / (2.f * h);

        const scalar evaluated = rq.derive_beta2();

        EXPECT_NEAR(numerical, evaluated, numerical * 0.01);
    }
}

// Test class for the derivative of bethe stopping power
class DerivativeOfBetheEquationValidation
    : public ::testing::TestWithParam<std::tuple<material<scalar>>> {};

TEST_P(DerivativeOfBetheEquationValidation, derivative_of_stopping_power) {

    // Interaction object
    interaction<scalar> Interactor;

    // Material
    const auto mat = std::get<0>(GetParam());

    // muon
    constexpr int pdg = pdg_particle::eMuon;

    // mass
    constexpr scalar m{105.7f * unit<scalar>::MeV};

    // charge
    const scalar q = -1.f;

    // Displacement for numerical differentiaion
    const scalar h = 1e-3f;

    // Iterate from 2 GeV to 100 GeV
    for (unsigned int i = 2u; i < 100; i++) {
        const scalar p = static_cast<scalar>(i) * detray::unit<scalar>::GeV;
        const scalar qop = q / p;

        const detail::relativistic_quantities<scalar> rq(m, qop, q);
        const detail::relativistic_quantities<scalar> rq1(m, qop + h, q);
        const detail::relativistic_quantities<scalar> rq2(m, qop - h, q);

        // Log term
        const scalar log_term1 = rq1.compute_bethe_bloch_log_term(mat);
        const scalar log_term2 = rq2.compute_bethe_bloch_log_term(mat);

        const scalar numerical_log_term = (log_term1 - log_term2) / (2.f * h);
        const scalar evaluated_log_term = rq.derive_bethe_bloch_log_term();

        EXPECT_NEAR(numerical_log_term, evaluated_log_term,
                    numerical_log_term * 0.01f);

        // delta half
        const scalar dhalf1 = rq1.compute_delta_half(mat);
        const scalar dhalf2 = rq2.compute_delta_half(mat);

        const scalar numerical_dhalf = (dhalf1 - dhalf2) / (2.f * h);
        const scalar evaluated_dhalf = rq.derive_delta_half(mat);

        EXPECT_NEAR(numerical_dhalf, evaluated_dhalf, numerical_dhalf * 0.01f);

        // Bethe equation
        const scalar bethe = Interactor.compute_bethe_bloch(mat, pdg, rq);
        const scalar bethe1 = Interactor.compute_bethe_bloch(mat, pdg, rq1);
        const scalar bethe2 = Interactor.compute_bethe_bloch(mat, pdg, rq2);

        const scalar numerical_bethe = (bethe1 - bethe2) / (2.f * h);

        const scalar evaluated_bethe =
            Interactor.derive_bethe_bloch(mat, pdg, rq, bethe);

        EXPECT_NEAR(numerical_bethe, evaluated_bethe, numerical_bethe * 0.01f);
    }
}

INSTANTIATE_TEST_SUITE_P(
    detray_material_BetheEquationDerivative,
    DerivativeOfBetheEquationValidation,
    ::testing::Values(hydrogen_gas<scalar>(), helium_gas<scalar>(),
                      isobutane<scalar>(), aluminium<scalar>(),
                      silicon<scalar>(), tungsten<scalar>(), gold<scalar>(),
                      cesium_iodide<scalar>(), silicon_with_ded<scalar>(),
                      cesium_iodide_with_ded<scalar>()));
