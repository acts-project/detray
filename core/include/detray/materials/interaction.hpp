/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/detail/relativistic_quantities.hpp"

namespace detray {

template <typename scalar_t>
struct interaction {

    using scalar_type = scalar_t;
    using relativistic_quantities =
        detail::relativistic_quantities<scalar_type>;

    DETRAY_HOST_DEVICE scalar_type compute_stopping_power(
        const detray::material<scalar_type>& mat, const int /*pdg*/,
        const scalar_type m, const scalar_type qOverP,
        const scalar_type q) const {

        const relativistic_quantities rq(m, qOverP, q);
        const scalar_t Ne{mat.molar_electron_density()};
        const scalar_t eps_per_length{rq.compute_epsilon_per_length(Ne)};
        if (eps_per_length <= 0.f) {
            return 0.f;
        }

        const scalar_t I{mat.mean_excitation_energy()};
        const scalar_t dhalf{rq.compute_delta_half(mat)};
        const scalar_t u{rq.compute_mass_term(constant<scalar_t>::m_e)};
        const scalar_t wmax{rq.compute_WMax(m)};
        // uses RPP2023 eq. 34.5 scaled from mass stopping power to linear
        // stopping power and multiplied with the material thickness to get a
        // total energy loss instead of an energy loss per length. the required
        // modification only change the prefactor which becomes identical to the
        // prefactor epsilon for the most probable value.
        const scalar_t running{math_ns::log(u * wmax / (I * I)) -
                               2.f * (rq.m_beta2 + dhalf)};
        return eps_per_length * running;
    }

    DETRAY_HOST_DEVICE scalar_type derive_stopping_power(
        const detray::material<scalar_type>& mat, const int /*pdg*/,
        const scalar_type m, const scalar_type qOverP, const scalar_type q) {

        const relativistic_quantities rq(m, qOverP, q);
        const scalar_t Ne{mat.molar_electron_density()};
        const scalar_t eps_per_length{rq.compute_epsilon_per_length(Ne)};
        if (eps_per_length <= 0.f) {
            return 0.f;
        }

        const scalar_t I{mat.mean_excitation_energy()};
        const scalar_t dhalf{rq.compute_delta_half(mat)};
        const scalar_t u{rq.compute_mass_term(constant<scalar_t>::m_e)};
        const scalar_t wmax{rq.compute_WMax(m)};
        // uses RPP2023 eq. 34.5 scaled from mass stopping power to linear
        // stopping power and multiplied with the material thickness to get a
        // total energy loss instead of an energy loss per length. the required
        // modification only change the prefactor which becomes identical to the
        // prefactor epsilon for the most probable value.
        const scalar_t running{math_ns::log(u * wmax / (I * I)) -
                               2.f * (rq.m_beta2 + dhalf)};

        // Get d(Beta)/d(qop)
        const scalar_t dBeta_dQop = rq.derive_beta();

        // W Max denominator
        const scalar_t denom = rq.compute_WMax_denominator(m);

        return -2.f * eps_per_length / rq.m_beta * dBeta_dQop *
               (running - (1.f + rq.m_betaGamma * rq.m_betaGamma) +
                2.f * rq.m_beta2 +
                constant<scalar_t>::m_e / m * rq.m_gamma * rq.m_betaGamma *
                    rq.m_betaGamma / denom);
    }

    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_bethe(
        const scalar_type path_segment,
        const detray::material<scalar_type>& mat, const int pdg,
        const scalar_type m, const scalar_type qOverP,
        const scalar_type q) const {

        return path_segment * compute_stopping_power(mat, pdg, m, qOverP, q);
    }

    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau(
        const scalar_type path_segment, const material<scalar_type>& mat,
        const int /*pdg*/, const scalar_type m, const scalar_type qOverP,
        const scalar_type q) const {

        const scalar_t I{mat.mean_excitation_energy()};
        const scalar_t Ne{mat.molar_electron_density()};
        const relativistic_quantities rq(m, qOverP, q);
        const scalar_t eps{rq.compute_epsilon(Ne, path_segment)};

        if (eps <= 0.f) {
            return 0.f;
        }

        const scalar_t dhalf{rq.compute_delta_half(mat)};
        const scalar_t t{rq.compute_mass_term(constant<scalar_t>::m_e)};
        // uses RPP2018 eq. 33.11
        const scalar_t running{math_ns::log(t / I) + math_ns::log(eps / I) +
                               0.2f - rq.m_beta2 - 2.f * dhalf};
        return eps * running;
    }

    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_fwhm(
        const scalar_type path_segment, const material<scalar_type>& mat,
        const int /*pdg*/, const scalar_type m, const scalar_type qOverP,
        const scalar_type q) const {
        const auto Ne = mat.molar_electron_density();
        const relativistic_quantities rq(m, qOverP, q);

        // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
        return 4.f * rq.compute_epsilon(Ne, path_segment);
    }

    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_sigma(
        const scalar_type path_segment, const material<scalar_type>& mat,
        const int pdg, const scalar_type m, const scalar_type qOverP,
        const scalar_type q) const {

        const scalar_t fwhm{compute_energy_loss_landau_fwhm(path_segment, mat,
                                                            pdg, m, qOverP, q)};

        return convert_landau_fwhm_to_gaussian_sigma(fwhm);
    }

    template <typename material_t>
    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_sigma_QOverP(
        const scalar_type path_segment, const material_t& mat, const int pdg,
        const scalar_type m, const scalar_type qOverP,
        const scalar_type q) const {

        const scalar_t sigmaE{compute_energy_loss_landau_sigma(
            path_segment, mat, pdg, m, qOverP, q)};

        //  var(q/p) = (d(q/p)/dE)² * var(E)
        // d(q/p)/dE = d/dE (q/sqrt(E²-m²))
        //           = q * -(1/2) * 1/p³ * 2E
        //  var(q/p) = (q/p)^4 * (q/beta)² * (1/q)^4 * var(E)
        //           = -q/p² E/p = -(q/p)² * 1/(q*beta) = -(q/p)² * (q/beta) /
        //           q²
        //           = (1/p)^4 * (q/beta)² * var(E)
        // do not need to care about the sign since it is only used squared
        const scalar_t pInv{qOverP / q};

        const relativistic_quantities rq(m, qOverP, q);

        return std::sqrt(rq.m_q2OverBeta2) * pInv * pInv * sigmaE;
    }

    DETRAY_HOST_DEVICE scalar_type compute_multiple_scattering_theta0(
        const scalar_type xOverX0, const int pdg, const scalar_type m,
        const scalar_type qOverP, const scalar_type q) const {

        // 1/p = q/(pq) = (q/p)/q
        const scalar_type momentumInv{std::abs(qOverP / q)};
        // q²/beta²; a smart compiler should be able to remove the unused
        // computations
        const relativistic_quantities rq(m, qOverP, q);
        const scalar_type q2OverBeta2{rq.m_q2OverBeta2};

        // if electron or positron
        if ((pdg == pdg_particle::eElectron) or
            (pdg == pdg_particle::ePositron)) {
            //@todo (Beomki): Not sure if we need this function. At least we
            // need to find the reference for this equation
            return theta0RossiGreisen(xOverX0, momentumInv, q2OverBeta2);
        } else {
            return theta0Highland(xOverX0, momentumInv, q2OverBeta2);
        }
    }

    private:
    /// Multiple scattering (mainly due to Coulomb interaction) for charged
    /// particles
    /// Original source: G. R. Lynch and O. I. Dahl, NIM.B58, 6
    DETRAY_HOST_DEVICE scalar_type
    theta0Highland(const scalar_type xOverX0, const scalar_type momentumInv,
                   const scalar_type q2OverBeta2) const {
        if (xOverX0 <= 0.f) {
            return 0.f;
        }

        // RPP2018 eq. 33.15 (treats beta and q² consistenly)
        const scalar_type t{std::sqrt(xOverX0 * q2OverBeta2)};
        // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
        //                          = 2 * log(sqrt(x/X0) * (q/beta))
        return 13.6f * unit<scalar_type>::MeV * momentumInv * t *
               (1.0f + 0.038f * 2.f * math_ns::log(t));
    }

    /// Multiple scattering theta0 for electrons.
    DETRAY_HOST_DEVICE scalar_type
    theta0RossiGreisen(const scalar_type xOverX0, const scalar_type momentumInv,
                       const scalar_type q2OverBeta2) const {
        if (xOverX0 <= 0.f) {
            return 0.f;
        }

        // TODO add source paper/ resource
        const scalar_type t{std::sqrt(xOverX0 * q2OverBeta2)};
        return 17.5f * unit<scalar_type>::MeV * momentumInv * t *
               (1.0f + 0.125f * math_ns::log10(10.0f * xOverX0));
    }

    /// Convert Landau full-width-half-maximum to an equivalent Gaussian sigma,
    ///
    /// Full-width-half-maximum for a Gaussian is given as
    ///
    ///     fwhm = 2 * sqrt(2 * log(2)) * sigma
    /// -> sigma = fwhm / (2 * sqrt(2 * log(2)))
    ///
    /// @todo: Add a unit test for this function
    DETRAY_HOST_DEVICE scalar_type
    convert_landau_fwhm_to_gaussian_sigma(const scalar_type fwhm) const {
        return 0.5f * constant<scalar_type>::inv_sqrt2 * fwhm /
               std::sqrt(constant<scalar_type>::ln2);
    }
};

}  // namespace detray
