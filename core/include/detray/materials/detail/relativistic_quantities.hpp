/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// Detray include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"

namespace detray::detail {

template <typename scalar_t>
struct relativistic_quantities {

    using scalar_type = scalar_t;

    // values from RPP2018 table 33.1
    // electron mass
    const scalar_type Me = 0.5109989461 * unit_constants::MeV;
    // Bethe formular prefactor. 1/mol unit is just a factor 1 here.
    const scalar_type K = 0.307075 * unit_constants::MeV * unit_constants::cm2;
    // Energy scale for plasma energy.
    const scalar_type PlasmaEnergyScale = 28.816 * unit_constants::eV;

    scalar_type m_q2OverBeta2 = 0.0;
    scalar_type m_beta2 = 0.0;
    scalar_type m_betaGamma = 0.0;
    scalar_type m_gamma = 0.0;

    relativistic_quantities(const scalar_type mass, const scalar_type qOverP,
                            const scalar_type q) {
        // beta²/q² = (p/E)²/q² = p²/(q²m² + q²p²) = 1/(q² + (m²(q/p)²)
        // q²/beta² = q² + m²(q/p)²
        m_q2OverBeta2 = q * q + (mass * qOverP) * (mass * qOverP);
        // 1/p = q/(qp) = (q/p)/q
        const auto mOverP = mass * std::abs(qOverP / q);
        const auto pOverM = scalar_type(1.0) / mOverP;
        // beta² = p²/E² = p²/(m² + p²) = 1/(1 + (m/p)²)
        m_beta2 = scalar_type(1.0) / (scalar_type(1.0) + mOverP * mOverP);
        // beta*gamma = (p/sqrt(m² + p²))*(sqrt(m² + p²)/m) = p/m
        m_betaGamma = pOverM;
        // gamma = sqrt(m² + p²)/m = sqrt(1 + (p/m)²)
        m_gamma = std::sqrt(scalar_type(1.0) + pOverM * pOverM);
    }

    /// Compute q/p derivative of beta².
    DETRAY_HOST_DEVICE inline scalar_type derive_beta2(
        const scalar_type qOverP) const {
        return -2 / (qOverP * m_gamma * m_gamma);
    }

    /// Compute the 2 * mass * (beta * gamma)² mass term.
    DETRAY_HOST_DEVICE inline scalar_type compute_mass_term(
        const scalar_type mass) const {
        return 2 * mass * m_betaGamma * m_betaGamma;
    }

    /// Compute mass term logarithmic derivative w/ respect to q/p.
    DETRAY_HOST_DEVICE inline scalar_type log_derive_mass_term(
        const scalar_type qOverP) const {
        // only need to compute d((beta*gamma)²)/(beta*gamma)²; rest cancels.
        return -2 / qOverP;
    }

    /// Compute the maximum energy transfer in a single collision.
    ///
    /// Uses RPP2018 eq. 33.4.
    DETRAY_HOST_DEVICE inline scalar_type compute_WMax(
        const scalar_type mass) const {
        const auto mfrac = Me / mass;
        const auto nominator = 2 * Me * m_betaGamma * m_betaGamma;
        const auto denonimator = 1.0 + 2 * m_gamma * mfrac + mfrac * mfrac;
        return nominator / denonimator;
    }

    /// Compute WMax logarithmic derivative w/ respect to q/p.
    DETRAY_HOST_DEVICE inline scalar_type log_derive_WMax(
        const scalar_type mass, const scalar_type qOverP) const {
        // this is (q/p) * (beta/q).
        // both quantities have the same sign and the product must always be
        // positive. we can thus reuse the known (unsigned) quantity (q/beta)².
        const auto a = std::abs(qOverP / std::sqrt(m_q2OverBeta2));
        // (m² + me²) / me = me (1 + (m/me)²)
        const auto b = Me * (1.0f + (mass / Me) * (mass / Me));
        return -2 * (a * b - 2 + m_beta2) / (qOverP * (a * b + 2));
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
        return 0.5 * K * molarElectronDensity * thickness * m_q2OverBeta2;
    }

    /// Compute epsilon logarithmic derivative w/ respect to q/p.
    DETRAY_HOST_DEVICE inline scalar_type log_derive_epsilon(
        const scalar_type qOverP) const {
        // only need to compute d(q²/beta²)/(q²/beta²); everything else cancels.
        return 2 / (qOverP * m_gamma * m_gamma);
    }

    /// Compute the density correction factor delta/2.
    ///
    /// Uses RPP2018 eq. 33.6 which is only valid for high energies.
    ///
    /// @todo Should we use RPP2018 eq. 33.7 instead w/ tabulated constants?
    DETRAY_HOST_DEVICE inline scalar_type compute_delta_half(
        const scalar_type meanExitationPotential,
        const scalar_type molarElectronDensity) const {
        // only relevant for very high ernergies; use arbitrary cutoff
        if (m_betaGamma < 10.0f) {
            return 0.0f;
        }
        // pre-factor according to RPP2019 table 33.1
        const auto plasmaEnergy =
            PlasmaEnergyScale * std::sqrt(molarElectronDensity);
        return std::log(m_betaGamma) +
               std::log(plasmaEnergy / meanExitationPotential) - 0.5f;
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