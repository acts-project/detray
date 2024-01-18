/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <cassert>

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

    scalar_type m_qOverP{0.f};
    scalar_type m_q2OverBeta2{0.f};
    scalar_type m_beta{0.f};
    scalar_type m_beta2{0.f};
    scalar_type m_betaGamma{0.f};
    scalar_type m_gamma{0.f};

    DETRAY_HOST_DEVICE
    relativistic_quantities(const scalar_type mass, const scalar_type qOverP,
                            const scalar_type q) {
        assert(qOverP != 0.f);
        assert(mass != 0.f);

        m_qOverP = qOverP;
        // beta²/q² = (p/E)²/q² = p²/(q²m² + q²p²) = 1/(q² + (m²(q/p)²)
        // q²/beta² = q² + m²(q/p)²
        m_q2OverBeta2 = q * q + (mass * qOverP) * (mass * qOverP);
        // 1/p = q/(qp) = (q/p)/q
        const scalar_type mOverP{
            mass * ((q == 0.f) ? std::abs(qOverP / q) : std::abs(qOverP))};
        const scalar_type pOverM{1.f / mOverP};
        // beta² = p²/E² = p²/(m² + p²) = 1/(1 + (m/p)²)
        m_beta2 = 1.f / (1.f + mOverP * mOverP);
        m_beta = math_ns::sqrt(m_beta2);
        // beta*gamma = (p/sqrt(m² + p²))*(sqrt(m² + p²)/m) = p/m
        m_betaGamma = pOverM;
        // gamma = sqrt(m² + p²)/m = sqrt(1 + (p/m)²)
        m_gamma = std::sqrt(1.f + pOverM * pOverM);
    }

    /// Compute the 2 * mass * (beta * gamma)² mass term.
    DETRAY_HOST_DEVICE inline constexpr scalar_type compute_mass_term(
        const scalar_type mass) const {
        return 2.f * mass * m_betaGamma * m_betaGamma;
    }

    /// Compute the maximum energy transfer in a single collision.
    ///
    /// Uses RPP2023 eq. 34.4.
    DETRAY_HOST_DEVICE inline constexpr scalar_type compute_WMax_denominator(
        const scalar_type mass) const {
        assert(mass != 0.f);
        const scalar_type mfrac{constant<scalar_type>::m_e / mass};
        return 1.f + 2.f * m_gamma * mfrac + mfrac * mfrac;
    }

    /// Compute the d(WMax_denominator)/d(qop)
    DETRAY_HOST_DEVICE inline constexpr scalar_type derive_WMax_denominator(
        const scalar_type mass) const {
        assert(mass != 0.f);
        assert(mass != 0.f);
        const scalar_type mfrac{constant<scalar_type>::m_e / mass};

        return 2.f * mfrac * m_betaGamma * m_gamma * m_gamma *
               this->derive_beta();
    }

    /// Compute the maximum energy transfer in a single collision.
    ///
    /// Uses RPP2023 eq. 34.4.
    DETRAY_HOST_DEVICE inline constexpr scalar_type compute_WMax(
        const scalar_type mass) const {
        const scalar_type nominator{2.f * constant<scalar_type>::m_e *
                                    m_betaGamma * m_betaGamma};
        return nominator / compute_WMax_denominator(mass);
    }

    // Compute d(WMax)/d(Qop)
    DETRAY_HOST_DEVICE inline constexpr scalar_type derive_WMax(
        const scalar_type mass) const {
        const scalar_type nominator{2.f * constant<scalar_type>::m_e *
                                    m_betaGamma * m_betaGamma};
        return nominator / compute_WMax_denominator(mass);
    }

    /// Compute epsilon per length where epsilon is defined at RPP2023
    /// eq. 34.12. Defined as
    ///
    ///     (K/2) * (Z/A) * rho * (q²/beta²)
    ///
    /// where (Z/A)*rho is the electron density in the material.
    DETRAY_HOST_DEVICE inline constexpr scalar_type compute_epsilon_per_length(
        const scalar_type molarElectronDensity) const {
        return 0.5f * K * molarElectronDensity * m_q2OverBeta2;
    }

    /// Compute epsilon energy pre-factor for RPP2023 eq. 34.12.
    ///
    /// Defined as
    ///
    ///     (K/2) * (Z/A) * rho * x * (q²/beta²)
    ///
    /// where (Z/A)*rho is the electron density in the material and x is the
    /// traversed length (thickness) of the material.
    DETRAY_HOST_DEVICE inline constexpr scalar_type compute_epsilon(
        const scalar_type molarElectronDensity,
        const scalar_type thickness) const {
        return compute_epsilon_per_length(molarElectronDensity) * thickness;
    }

    /// Compute the density correction factor delta/2.
    DETRAY_HOST_DEVICE inline scalar_type compute_delta_half(
        const material<scalar_type>& mat) const {

        if (!mat.has_density_effect_data()) {
            // Apply the approximated density effector correction to tracks with
            // beta*gamma > 10
            //
            // @NOTE The approximated function is only accurate for beta*gamma >
            // 100 therefore, the cutoff with 10 might not be a good choice. The
            // cutoff also introduces a step where the energy loss function is
            // not smooth anymore.
            //
            // For further discussion, please follow the ATLAS JIRA Tickets:
            // ATLASRECTS-3144 and ATLASRECTS-7586 (ATLAS-restricted)
            if (m_betaGamma < 10.f) {
                return 0.f;
            }
            // Equation 34.6 of PDG2022
            // @NOTE A factor of 1000 is required to convert the unit of density
            // (mm^-3 to cm^-3)
            const scalar_type plasmaEnergy{
                PlasmaEnergyScale *
                std::sqrt(1000.f * mat.molar_electron_density())};
            return math_ns::log(m_betaGamma * plasmaEnergy /
                                mat.mean_excitation_energy()) -
                   0.5f;
        } else {
            const auto& density = mat.density_effect_data();

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

    /// Derive the density correction factor delta/2.
    DETRAY_HOST_DEVICE inline scalar_type derive_delta_half(
        const material<scalar_type>& mat) const {

        if (!mat.has_density_effect_data()) {
            return this->derive_betaGamma() / m_betaGamma;
        } else {
            const auto& density = mat.density_effect_data();

            // const scalar_type cden{density.get_C_density()};
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
                delta = 2.f / m_betaGamma * this->derive_betaGamma();
                if (x < x1den)
                    delta +=
                        aden * mden * math_ns::pow(x1den - x, mden - 1) *
                        (-1.f / (m_betaGamma * constant<scalar_type>::ln10)) *
                        this->derive_betaGamma();
            }

            return 0.5f * delta;
        }
    }

    /// Compute derivative of beta w.r.t q/p
    ///
    /// Used the relation of [p = gamma * m * v]:
    /// -> dp/d(beta) = mc*gamma*(1+beta^2*gamma^2)
    ///
    /// d(beta)/d(qop) = -beta/[qop * (1+beta^2*gamma^2)]
    DETRAY_HOST_DEVICE inline constexpr scalar_type derive_beta() const {
        assert(m_qOverP != 0.f);
        assert(m_betaGamma != 0.f);
        return -1.f * m_beta / (m_qOverP * (1.f + m_betaGamma * m_betaGamma));
    }

    /// Compute derivative of betagamma w.r.t. q/p
    ///
    /// d(betagamma)/dqop = d(beta)/dqop * gamma + d(gamma)/dqop * beta
    ///
    /// Since d(gamma)/dqop = d(beta)/dqop * gamma * (beta / 1 - beta^2)
    ///
    /// d(betagamma)/dqop = gamma^3 * dbeta/dqop
    DETRAY_HOST_DEVICE inline constexpr scalar_type derive_betaGamma() const {
        assert(m_qOverP != 0.f);
        assert(m_betaGamma != 0.f);
        return m_gamma * m_gamma * m_gamma * this->derive_beta();
    }
};

}  // namespace detray::detail
