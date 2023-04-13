/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/material.hpp"

namespace detray::detail {

template <typename scalar_t>
struct relativistic_quantities {

    using scalar_type = scalar_t;

    // values from RPP2018 table 33.1
    // Bethe formular prefactor. 1/mol unit is just a factor 1 here.
    static constexpr scalar_type K{
        0.307075f * (unit<scalar_type>::MeV * unit<scalar_type>::cm2)};
    // Energy scale for plasma energy.
    static constexpr scalar_type PlasmaEnergyScale{28.816f *
                                                   unit<scalar_type>::eV};

    scalar_type m_q2OverBeta2{0.f};
    scalar_type m_beta2{0.f};
    scalar_type m_betaGamma{0.f};
    scalar_type m_gamma{0.f};

    DETRAY_HOST_DEVICE
    relativistic_quantities(const scalar_type mass, const scalar_type qOverP,
                            const scalar_type q) {
        // beta²/q² = (p/E)²/q² = p²/(q²m² + q²p²) = 1/(q² + (m²(q/p)²)
        // q²/beta² = q² + m²(q/p)²
        m_q2OverBeta2 = q * q + (mass * qOverP) * (mass * qOverP);
        // 1/p = q/(qp) = (q/p)/q
        const scalar_type mOverP{mass * std::abs(qOverP / q)};
        const scalar_type pOverM{1.f / mOverP};
        // beta² = p²/E² = p²/(m² + p²) = 1/(1 + (m/p)²)
        m_beta2 = 1.f / (1.f + mOverP * mOverP);
        // beta*gamma = (p/sqrt(m² + p²))*(sqrt(m² + p²)/m) = p/m
        m_betaGamma = pOverM;
        // gamma = sqrt(m² + p²)/m = sqrt(1 + (p/m)²)
        m_gamma = std::sqrt(1.f + pOverM * pOverM);
    }

    /// Compute q/p derivative of beta².
    DETRAY_HOST_DEVICE inline scalar_type derive_beta2(
        const scalar_type qOverP) const {
        return -2.f / (qOverP * m_gamma * m_gamma);
    }

    /// Compute the 2 * mass * (beta * gamma)² mass term.
    DETRAY_HOST_DEVICE inline scalar_type compute_mass_term(
        const scalar_type mass) const {
        return 2.f * mass * m_betaGamma * m_betaGamma;
    }

    /// Compute mass term logarithmic derivative w/ respect to q/p.
    DETRAY_HOST_DEVICE inline scalar_type log_derive_mass_term(
        const scalar_type qOverP) const {
        // only need to compute d((beta*gamma)²)/(beta*gamma)²; rest cancels.
        return -2.f / qOverP;
    }

    /// Compute the maximum energy transfer in a single collision.
    ///
    /// Uses RPP2018 eq. 33.4.
    DETRAY_HOST_DEVICE inline scalar_type compute_WMax(
        const scalar_type mass) const {
        const scalar_type mfrac{constant<scalar_type>::m_e / mass};
        const scalar_type nominator{2.f * constant<scalar_type>::m_e *
                                    m_betaGamma * m_betaGamma};
        const scalar_type denonimator{1.f + 2.f * m_gamma * mfrac +
                                      mfrac * mfrac};
        return nominator / denonimator;
    }

    /// Compute WMax logarithmic derivative w/ respect to q/p.
    DETRAY_HOST_DEVICE inline scalar_type log_derive_WMax(
        const scalar_type mass, const scalar_type qOverP) const {
        // this is (q/p) * (beta/q).
        // both quantities have the same sign and the product must always be
        // positive. we can thus reuse the known (unsigned) quantity (q/beta)².
        const scalar_type a{std::abs(qOverP / std::sqrt(m_q2OverBeta2))};
        // (m² + me²) / me = me (1 + (m/me)²)
        const scalar_type rel_mass{mass / constant<scalar_type>::m_e};
        const scalar_type b =
            constant<scalar_type>::m_e * (1.0f + rel_mass * rel_mass);
        return -2.f * (a * b - 2.f + m_beta2) / (qOverP * (a * b + 2.f));
    }

    /// Compute epsilon energy pre-factor for RPP2018 eq. 33.11.
    ///
    /// Defined as
    ///
    ///     (K/2) * (Z/A)*rho * x * (q²/beta²)
    ///
    /// where (Z/A)*rho is the electron density in the material and x is the
    /// traversed length (thickness) of the material.
    DETRAY_HOST_DEVICE inline scalar_type compute_epsilon(
        const scalar_type molarElectronDensity,
        const scalar_type thickness) const {
        return 0.5f * K * molarElectronDensity * thickness * m_q2OverBeta2;
    }

    /// Compute epsilon logarithmic derivative w/ respect to q/p.
    DETRAY_HOST_DEVICE inline scalar_type log_derive_epsilon(
        const scalar_type qOverP) const {
        // only need to compute d(q²/beta²)/(q²/beta²); everything else cancels.
        return 2.f / (qOverP * m_gamma * m_gamma);
    }

    /// Compute the density correction factor delta/2.
    DETRAY_HOST_DEVICE inline scalar_type compute_delta_half(
        const material<scalar_type>& mat) const {

        const auto& density = mat.density_effect_data();

        // When the denstiy effect data is not provided (RPP2018 eq. 33.6 only
        // valid for high energies)
        if (density == density_effect_data<scalar_type>{}) {
            const scalar_type mean_exitation_energy =
                mat.mean_excitation_energy();
            const scalar_type molar_electron_density =
                mat.molar_electron_density();

            // only relevant for very high ernergies; use arbitrary cutoff
            if (m_betaGamma < 10.f) {
                return 0.f;
            }
            // pre-factor according to RPP2019 table 33.1
            const scalar_type plasmaEnergy{
                PlasmaEnergyScale * std::sqrt(1000.f * molar_electron_density)};
            return math_ns::log(m_betaGamma) +
                   math_ns::log(plasmaEnergy / mean_exitation_energy) - 0.5f;
        }
        // When the denstiy effect data is provided (RPP2018 eq. 33.7)
        else {
            const scalar_type cden{density.get_C_density()};
            const scalar_type mden{density.get_M_density()};
            const scalar_type aden{density.get_A_density()};
            const scalar_type x0den{density.get_X0_density()};
            const scalar_type x1den{density.get_X1_density()};

            const scalar_type x{math_ns::log10(m_betaGamma)};

            scalar_type delta;

            // From Geant4
            // processes/electromagnetic/lowenergy/src/G4hBetheBlochModel.cc
            if (x < x0den) {
                delta = 0.f;

            } else {
                delta = 2.f * constant<scalar_type>::ln10 * x - cden;
                if (x < x1den)
                    delta += aden * math_ns::pow((x1den - x), mden);
            }

            return 0.5f * delta;
        }
    }

    /// Compute derivative w/ respect to q/p for the density correction.
    DETRAY_HOST_DEVICE inline scalar_type derive_delta_half(
        const scalar_type qOverP) const {
        // original equation is of the form
        //     log(beta*gamma) + log(eplasma/I) - 1/2
        // which the resulting derivative as
        //     d(beta*gamma)/(beta*gamma)
        return (m_betaGamma < 10.0f) ? 0.0f : (-1.0f / qOverP);
    }
};

}  // namespace detray::detail