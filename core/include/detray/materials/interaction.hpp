/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/detail/relativistic_quantities.hpp"

namespace detray {

template <typename scalar_t>
struct interaction {

    using scalar_type = scalar_t;
    using relativistic_quantities =
        detail::relativistic_quantities<scalar_type>;

    // @returns the total stopping power
    DETRAY_HOST_DEVICE scalar_type
    compute_stopping_power(const detray::material<scalar_type>& mat,
                           const pdg_particle<scalar_type>& ptc,
                           const relativistic_quantities& rq) const {

        scalar_type stopping_power{0.f};
        stopping_power += compute_bethe_bloch(mat, ptc, rq);
        // stopping_power += compute_bremsstralung(mat, ptc, rq);

        return stopping_power;
    }

    DETRAY_HOST_DEVICE scalar_type
    compute_bethe_bloch(const detray::material<scalar_type>& mat,
                        const pdg_particle<scalar_type>& /*ptc*/,
                        const relativistic_quantities& rq) const {

        const scalar_type eps_per_length{rq.compute_epsilon_per_length(mat)};
        if (eps_per_length <= 0.f) {
            return 0.f;
        }

        const scalar_type dhalf{rq.compute_delta_half(mat)};
        const scalar_type A = rq.compute_bethe_bloch_log_term(mat);
        const scalar_type running{A - rq.m_beta2 - dhalf};
        return 2.f * eps_per_length * running;
    }

    DETRAY_HOST_DEVICE scalar_type derive_bethe_bloch(
        const detray::material<scalar_type>& mat,
        const pdg_particle<scalar_type>& /*ptc*/,
        const relativistic_quantities& rq, const scalar_type bethe) const {

        // (K/2) * (Z/A) * z^2 / beta^2 * density
        const scalar_type eps_per_length{rq.compute_epsilon_per_length(mat)};
        if (eps_per_length <= 0.f) {
            return 0.f;
        }

        /*-----------------------------------------------------------------------
        * Calculation of d(-dE/dx)/dqop
        *
        * d(-dE/dx)/dqop = 2/(qop * gamma^2) * (-dE/dx)
        *                  + 2 * (eps/x) * [dA/dqop - dB/dqop - dC/dqop]
        *
        * where
        * A = 1/2 ln ( 2m_e c^2 beta^2 gamma^2 W_max/ I^2)
        * B = beta^2
        * C = delta/2
        *
        * dA/dqop = -1 / (2 * qop) * [4 - W_max/ (gamma M c^2) ]
        * dB/dqop = -2 * beta^2/(qop gamma^2)
        *
        * dC/dqop = 1/2*(-2/qop) if (x>x_1)
        *         = 1/2*(-2/qop + ak/(qop ln10)*(x_1 - x)^(k-1)) if (x_0<x<x_1)
        *         = 0 if (x<x_0) (for nonconductors)
        ------------------------------------------------------------------------*/

        const scalar_type first_term =
            2.f / (rq.m_qOverP * rq.m_gamma2) * bethe;

        const scalar_type dAdqop = rq.derive_bethe_bloch_log_term();
        const scalar_type dBdqop = rq.derive_beta2();
        const scalar_type dCdqop = rq.derive_delta_half(mat);

        const scalar_type second_term =
            2.f * eps_per_length * (dAdqop - dBdqop - dCdqop);

        return first_term + second_term;
    }

    DETRAY_HOST_DEVICE scalar_type
    compute_energy_loss_bethe_bloch(const scalar_type path_segment,
                                    const detray::material<scalar_type>& mat,
                                    const pdg_particle<scalar_type>& ptc,
                                    const relativistic_quantities& rq) const {

        return path_segment * compute_bethe_bloch(mat, ptc, rq);
    }

    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau(
        const scalar_type path_segment, const material<scalar_type>& mat,
        const pdg_particle<scalar_type>& /*ptc*/,
        const relativistic_quantities& rq) const {

        const scalar_type I{mat.mean_excitation_energy()};
        const scalar_type eps{rq.compute_epsilon(mat, path_segment)};

        if (eps <= 0.f) {
            return 0.f;
        }

        const scalar_type dhalf{rq.compute_delta_half(mat)};
        const scalar_type t{rq.compute_mass_term(constant<scalar_type>::m_e)};
        // uses RPP2018 eq. 33.11
        const scalar_type running{math::log(t / I) + math::log(eps / I) + 0.2f -
                                  rq.m_beta2 - 2.f * dhalf};
        return eps * running;
    }

    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_fwhm(
        const scalar_type path_segment, const material<scalar_type>& mat,
        const pdg_particle<scalar_type>& /*ptc*/,
        const relativistic_quantities& rq) const {

        // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
        return 4.f * rq.compute_epsilon(mat, path_segment);
    }

    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_sigma(
        const scalar_type path_segment, const material<scalar_type>& mat,
        const pdg_particle<scalar_type>& ptc,
        const relativistic_quantities& rq) const {

        const scalar_type fwhm{
            compute_energy_loss_landau_fwhm(path_segment, mat, ptc, rq)};

        return convert_landau_fwhm_to_gaussian_sigma(fwhm);
    }

    template <typename material_t>
    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_sigma_QOverP(
        const scalar_type path_segment, const material_t& mat,
        const pdg_particle<scalar_type>& ptc,
        const relativistic_quantities& rq) const {

        const scalar_type sigmaE{
            compute_energy_loss_landau_sigma(path_segment, mat, ptc, rq)};

        //  var(q/p) = (d(q/p)/dE)² * var(E)
        // d(q/p)/dE = d/dE (q/sqrt(E²-m²))
        //           = q * -(1/2) * 1/p³ * 2E
        //  var(q/p) = (q/p)^4 * (q/beta)² * (1/q)^4 * var(E)
        //           = -q/p² E/p = -(q/p)² * 1/(q*beta) = -(q/p)² * (q/beta)
        //           / q² = (1/p)^4 * (q/beta)² * var(E)
        // do not need to care about the sign since it is only used squared
        const scalar_type pInv{rq.m_qOverP / ptc.charge()};
        return math::sqrt(rq.m_q2OverBeta2) * pInv * pInv * sigmaE;
    }

    DETRAY_HOST_DEVICE scalar_type compute_multiple_scattering_theta0(
        const scalar_type xOverX0, const pdg_particle<scalar_type>& ptc,
        const relativistic_quantities& rq) const {

        // 1/p = q/(pq) = (q/p)/q
        const scalar_type momentumInv{math::fabs(rq.m_qOverP / ptc.charge())};
        // q²/beta²; a smart compiler should be able to remove the unused
        // computations

        // if electron or positron
        if ((ptc.pdg_num() == 11) or (ptc.pdg_num() == -11)) {
            //@todo (Beomki): Not sure if we need this function. At least we
            // need to find the reference for this equation
            return theta0RossiGreisen(xOverX0, momentumInv, rq.m_q2OverBeta2);
        } else {
            return theta0Highland(xOverX0, momentumInv, rq.m_q2OverBeta2);
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
        const scalar_type t{math::sqrt(xOverX0 * q2OverBeta2)};
        // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
        //                          = 2 * log(sqrt(x/X0) * (q/beta))
        return 13.6f * unit<scalar_type>::MeV * momentumInv * t *
               (1.0f + 0.038f * 2.f * math::log(t));
    }

    /// Multiple scattering theta0 for electrons.
    DETRAY_HOST_DEVICE scalar_type
    theta0RossiGreisen(const scalar_type xOverX0, const scalar_type momentumInv,
                       const scalar_type q2OverBeta2) const {
        if (xOverX0 <= 0.f) {
            return 0.f;
        }

        // TODO add source paper/ resource
        const scalar_type t{math::sqrt(xOverX0 * q2OverBeta2)};
        return 17.5f * unit<scalar_type>::MeV * momentumInv * t *
               (1.0f + 0.125f * math::log10(10.0f * xOverX0));
    }

    /// Convert Landau full-width-half-maximum to an equivalent Gaussian
    /// sigma,
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
               math::sqrt(constant<scalar_type>::ln2);
    }
};

}  // namespace detray
