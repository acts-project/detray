/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/materials/detail/relativistic_quantities.hpp"

namespace detray {

template <typename scalar_t>
struct interaction {

    using scalar_type = scalar_t;
    using matrix_operator = standard_matrix_operator<scalar_type>;
    using relativistic_quantities =
        detail::relativistic_quantities<scalar_type>;

    template <typename material_t>
    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_bethe(
        const line_plane_intersection& is, const material_t& mat,
        const int /*pdg*/, const scalar_type m, const scalar_type qOverP,
        const scalar_type q = unit_constants::e) const {

        // return early in case of vacuum or zero thickness
        if (not mat) {
            return scalar_type(0.);
        }

        const auto I = mat.get_material().mean_excitation_energy();
        const auto Ne = mat.get_material().molar_electron_density();
        const auto path_segment = mat.path_segment(is);
        const relativistic_quantities rq(m, qOverP, q);
        const auto eps = rq.compute_epsilon(Ne, path_segment);
        const auto dhalf = rq.compute_delta_half(mat.get_material());
        const auto u = rq.compute_mass_term(rq.Me);
        const auto wmax = rq.compute_WMax(m);
        // uses RPP2018 eq. 33.5 scaled from mass stopping power to linear
        // stopping power and multiplied with the material thickness to get a
        // total energy loss instead of an energy loss per length. the required
        // modification only change the prefactor which becomes identical to the
        // prefactor epsilon for the most probable value.
        const auto running =
            std::log(u / I) + std::log(wmax / I) - 2. * rq.m_beta2 - 2. * dhalf;
        return eps * running;
    }

    template <typename material_t>
    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau(
        const line_plane_intersection& is, const material_t& mat,
        const int /*pdg*/, const scalar_type m, const scalar_type qOverP,
        const scalar_type q = unit_constants::e) const {

        // return early in case of vacuum or zero thickness
        if (not mat) {
            return scalar_type(0.);
        }

        const auto I = mat.get_material().mean_excitation_energy();
        const auto Ne = mat.get_material().molar_electron_density();
        const auto path_segment = mat.path_segment(is);
        const relativistic_quantities rq(m, qOverP, q);
        const auto eps = rq.compute_epsilon(Ne, path_segment);
        const auto dhalf = rq.compute_delta_half(mat.get_material());
        const auto t = rq.compute_mass_term(rq.Me);
        // uses RPP2018 eq. 33.11
        const auto running = std::log(t / I) + std::log(eps / I) +
                             scalar_type(0.2) - rq.m_beta2 -
                             scalar_type(2.) * dhalf;
        return eps * running;
    }

    template <typename material_t>
    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_fwhm(
        const line_plane_intersection& is, const material_t& mat,
        const int /*pdg*/, const scalar_type m, const scalar_type qOverP,
        const scalar_type q = unit_constants::e) const {
        const auto Ne = mat.get_material().molar_electron_density();
        const auto path_segment = mat.path_segment(is);
        const relativistic_quantities rq(m, qOverP, q);

        // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
        return scalar_type(4.) * rq.compute_epsilon(Ne, path_segment);
    }

    template <typename material_t>
    DETRAY_HOST_DEVICE scalar_type compute_energy_loss_landau_sigma_QOverP(
        const line_plane_intersection& is, const material_t& mat, const int pdg,
        const scalar_type m, const scalar_type qOverP,
        const scalar_type q = unit_constants::e) const {

        // return early in case of vacuum or zero thickness
        if (not mat) {
            return scalar_type(0.);
        }
        /*
        const auto Ne = mat.get_material().molar_electron_density();
        const auto path_segment = mat.path_segment(is);
        const relativistic_quantities rq(m, qOverP, q);
        // the Landau-Vavilov fwhm is 4*eps (see RPP2018 fig. 33.7)
        const auto fwhm = 4 * rq.compute_epsilon(Ne, path_segment);
        */
        const relativistic_quantities rq(m, qOverP, q);
        const auto fwhm =
            compute_energy_loss_landau_fwhm(is, mat, pdg, m, qOverP, q);
        const auto sigmaE = convert_landau_fwhm_to_gaussian_sigma(fwhm);
        //  var(q/p) = (d(q/p)/dE)² * var(E)
        // d(q/p)/dE = d/dE (q/sqrt(E²-m²))
        //           = q * -(1/2) * 1/p³ * 2E
        //  var(q/p) = (q/p)^4 * (q/beta)² * (1/q)^4 * var(E)
        //           = -q/p² E/p = -(q/p)² * 1/(q*beta) = -(q/p)² * (q/beta) /
        //           q²
        //           = (1/p)^4 * (q/beta)² * var(E)
        // do not need to care about the sign since it is only used squared
        const auto pInv = qOverP / q;
        return std::sqrt(rq.m_q2OverBeta2) * pInv * pInv * sigmaE;
    }

    template <typename material_t>
    DETRAY_HOST_DEVICE scalar_type compute_multiple_scattering_theta0(
        const line_plane_intersection& is, const material_t& mat, const int pdg,
        const scalar_type m, const scalar_type qOverP,
        const scalar_type q = unit_constants::e) const {
        // return early in case of vacuum or zero thickness
        if (not mat) {
            return scalar_type(0.);
        }

        // relative radiation length
        const scalar_type xOverX0 = mat.path_segment_in_X0(is);
        // 1/p = q/(pq) = (q/p)/q
        const scalar_type momentumInv = std::abs(qOverP / q);
        // q²/beta²; a smart compiler should be able to remove the unused
        // computations
        const relativistic_quantities rq(m, qOverP, q);
        const scalar_type q2OverBeta2 = rq.m_q2OverBeta2;

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
    /// Multiple scattering theta0 for minimum ionizing particles.
    DETRAY_HOST_DEVICE scalar_type
    theta0Highland(const scalar_type xOverX0, const scalar_type momentumInv,
                   const scalar_type q2OverBeta2) const {
        // RPP2018 eq. 33.15 (treats beta and q² consistenly)
        const scalar_type t = std::sqrt(xOverX0 * q2OverBeta2);
        // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
        //                          = 2 * log(sqrt(x/X0) * (q/beta))
        return 13.6 * unit_constants::MeV * momentumInv * t *
               (1.0 + 0.038 * 2. * std::log(t));
    }

    /// Multiple scattering theta0 for electrons.
    DETRAY_HOST_DEVICE scalar_type
    theta0RossiGreisen(const scalar_type xOverX0, const scalar_type momentumInv,
                       const scalar_type q2OverBeta2) const {
        // TODO add source paper/ resource
        const scalar_type t = std::sqrt(xOverX0 * q2OverBeta2);
        return 17.5 * unit_constants::MeV * momentumInv * t *
               (1.0 + 0.125 * std::log10(10.0 * xOverX0));
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
        return fwhm / (2 * std::sqrt(2 * std::log(2.0)));
    }
};

}  // namespace detray