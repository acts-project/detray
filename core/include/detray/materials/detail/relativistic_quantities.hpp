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
#include "detray/materials/material.hpp"

namespace detray::detail {

template <typename scalar_t>
struct relativistic_quantities {

    using scalar_type = scalar_t;

    // values from RPP2018 table 33.1
    // electron mass
    const scalar_type Me = 0.5109989461 * unit<scalar_type>::MeV;
    // Bethe formular prefactor. 1/mol unit is just a factor 1 here.
    const scalar_type K =
        0.307075 * unit<scalar_type>::MeV * unit<scalar_type>::cm2;
    // Energy scale for plasma energy.
    const scalar_type PlasmaEnergyScale = 28.816 * unit<scalar_type>::eV;

    scalar_type m_q2OverBeta2 = 0.0;
    scalar_type m_beta2 = 0.0;
    scalar_type m_betaGamma = 0.0;
    scalar_type m_gamma = 0.0;

    DETRAY_HOST_DEVICE
    relativistic_quantities(const scalar_type mass, const scalar_type qOverP,
                            const scalar_type q) {
        // beta²/q² = (p/E)²/q² = p²/(q²m² + q²p²) = 1/(q² + (m²(q/p)²)
        // q²/beta² = q² + m²(q/p)²
        m_q2OverBeta2 = q * q + (mass * qOverP) * (mass * qOverP);
        // 1/p = q/(qp) = (q/p)/q
        const scalar_type mOverP = mass * std::abs(qOverP / q);
        const scalar_type pOverM = scalar_type(1.0) / mOverP;
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
    DETRAY_HOST_DEVICE inline scalar_type compute_delta_half(
        const material<scalar_type>& mat) const {

        const auto& density = mat.density_effect_data();

        // When the denstiy effect data is not provided (RPP2018 eq. 33.6 only
        // valid for high energies)
        if (density == density_effect_data<scalar_type>(0, 0, 0, 0, 0, 0, 0)) {
            const scalar_type mean_exitation_energy =
                mat.mean_excitation_energy();
            const scalar_type molar_electron_density =
                mat.molar_electron_density();

            // only relevant for very high ernergies; use arbitrary cutoff
            if (m_betaGamma < scalar_type(10.0)) {
                return scalar_type(0.);
            }
            // pre-factor according to RPP2019 table 33.1
            const scalar_type plasmaEnergy =
                PlasmaEnergyScale * std::sqrt(molar_electron_density);
            return std::log(m_betaGamma) +
                   std::log(plasmaEnergy / mean_exitation_energy) - 0.5f;
        }
        // When the denstiy effect data is provided (RPP2018 eq. 33.7)
        else {
            const scalar_type cden = density.get_C_density();
            const scalar_type mden = density.get_M_density();
            const scalar_type aden = density.get_A_density();
            const scalar_type x0den = density.get_X0_density();
            const scalar_type x1den = density.get_X1_density();

            const scalar_type x = std::log10(m_betaGamma);

            scalar_type delta;

            // From Geant4
            // processes/electromagnetic/lowenergy/src/G4hBetheBlochModel.cc
            if (x < x0den) {
                delta = 0.0;

            } else {
                delta = scalar_type(2.) * std::log(10.) * x - cden;
                if (x < x1den)
                    delta += aden * std::pow((x1den - x), mden);
            }

            return delta / scalar_type(2.);
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