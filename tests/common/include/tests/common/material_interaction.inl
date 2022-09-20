/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/actors/resetter.hpp"
#include "detray/propagator/actors/scattering_simulator.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/inspectors.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using transform3 = __plugin::transform3<scalar>;
using matrix_operator = typename transform3::matrix_actor;

// Test class for MUON energy loss with Bethe function
// Input tuple: < material / energy / expected output from
// https://pdg.lbl.gov/2022/AtomicNuclearProperties for Muon dEdX and range >
class EnergyLossBetheValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, scalar, scalar>> {};

// This tests the material functionalities
TEST_P(EnergyLossBetheValidation, bethe_energy_loss) {

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // H2 liquid with a unit thickness
    material_slab<scalar> slab(std::get<0>(GetParam()), 1 * unit_constants::cm);

    // muon
    const int pdg = pdg_particle::eMuon;

    // mass
    const scalar m = 105.7 * unit_constants::MeV;

    // qOverP
    const scalar qOverP = -1. / std::get<1>(GetParam());

    // Bethe Stopping power in MeV * cm^2 / g
    const scalar dEdx =
        I.compute_energy_loss_bethe(is, slab, pdg, m, qOverP, -1.) /
        slab.path_segment(is) / slab.get_material().mass_density() /
        (unit_constants::MeV * unit_constants::cm2 / unit_constants::g);

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - dEdx) / dEdx < 0.05);
}

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      0.1003 * unit_constants::GeV, 6.539)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      1.101 * unit_constants::GeV, 4.182)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      10.11 * unit_constants::GeV, 4.777)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(),
                                      100.1 * unit_constants::GeV, 5.305)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      0.1003 * unit_constants::GeV, 3.082)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      1.101 * unit_constants::GeV, 2.133)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      10.11 * unit_constants::GeV, 2.768)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(),
                                      100.1 * unit_constants::GeV, 3.188)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      0.1003 * unit_constants::GeV, 2.533)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      1.101 * unit_constants::GeV, 1.744)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      10.11 * unit_constants::GeV, 2.097)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(),
                                      100.1 * unit_constants::GeV, 2.360)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_0p1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      0.1003 * unit_constants::GeV, 2.608)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      1.101 * unit_constants::GeV, 1.803)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_10GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      10.11 * unit_constants::GeV, 2.177)));

INSTANTIATE_TEST_SUITE_P(
    Bethe_100GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      100.1 * unit_constants::GeV, 2.451)));

// Test class for MUON energy loss with Landau function
// Input tuple: < material / energy / expected energy loss  / expected fwhm  >
class EnergyLossLandauValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, scalar, scalar, scalar>> {};

TEST_P(EnergyLossLandauValidation, landau_energy_loss) {

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // H2 liquid with a unit thickness
    material_slab<scalar> slab(std::get<0>(GetParam()),
                               0.17 * unit_constants::cm);

    // muon
    const int pdg = pdg_particle::eMuon;

    // mass
    const scalar m = 105.7 * unit_constants::MeV;

    // qOverP
    const scalar qOverP = -1. / std::get<1>(GetParam());

    // Landau Energy loss in MeV
    const scalar dE =
        I.compute_energy_loss_landau(is, slab, pdg, m, qOverP, -1) /
        (unit_constants::MeV);

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - dE) / dE < 0.05);

    // Landau Energy loss Fluctuation
    const scalar fwhm =
        I.compute_energy_loss_landau_fwhm(is, slab, pdg, m, qOverP, -1) /
        (unit_constants::MeV);

    // Check if difference is within 10% error
    EXPECT_TRUE(std::abs(std::get<3>(GetParam()) - fwhm) / fwhm < 0.1);
}

// Expected output from Fig 33.7 in RPP2018
INSTANTIATE_TEST_SUITE_P(
    Landau_10GeV_Silicon, EnergyLossLandauValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(),
                                      10. * unit_constants::GeV, 0.525, 0.13)));

// Material interaction test with telescope Geometry
TEST(material_interaction, telescope_geometry_energy_loss) {

    vecmem::host_memory_resource host_mr;

    using ln_stepper_t = line_stepper<transform3>;

    typename ln_stepper_t::free_track_parameters_type default_trk{
        {0, 0, 0}, 0, {1, 0, 0}, -1};

    // Build from given module positions
    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};

    const auto mat = silicon_tml<scalar>();
    const scalar thickness = 0.17 * unit_constants::cm;

    const auto det = create_telescope_detector(
        host_mr, positions, ln_stepper_t(),
        typename ln_stepper_t::state{default_trk}, 20. * unit_constants::mm,
        20. * unit_constants::mm, mat, thickness);

    using navigator_t = navigator<decltype(det)>;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_default_policy;
    using stepper_t = line_stepper<transform3, constraints_t, policy_t>;
    using interactor_t = pointwise_material_interactor<transform3>;
    using actor_chain_t =
        actor_chain<dtuple, propagation::print_inspector, pathlimit_aborter,
                    parameter_transporter<transform3>, interactor_t,
                    resetter<transform3>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator is built from the stepper and navigator
    propagator_t p({}, {});

    const scalar q = -1.;
    const scalar iniP = 10 * unit_constants::GeV;

    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = 0.;
    getter::element(bound_vector, e_bound_loc1, 0) = 0.;
    getter::element(bound_vector, e_bound_phi, 0) = 0.;
    getter::element(bound_vector, e_bound_theta, 0) = M_PI / 2.;
    getter::element(bound_vector, e_bound_qoverp, 0) = q / iniP;
    getter::element(bound_vector, e_bound_time, 0) = 0.;
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    // bound track parameter
    const bound_track_parameters<transform3> bound_param(0, bound_vector,
                                                         bound_cov);

    propagation::print_inspector::state print_insp_state{};
    pathlimit_aborter::state aborter_state{};
    parameter_transporter<transform3>::state bound_updater{};
    interactor_t::state interactor_state{};
    resetter<transform3>::state resetter_state{};

    // Create actor states tuples
    actor_chain_t::state actor_states =
        std::tie(print_insp_state, aborter_state, bound_updater,
                 interactor_state, resetter_state);

    propagator_t::state state(bound_param, det, actor_states);

    // Propagate the entire detector
    ASSERT_TRUE(p.propagate(state))
        << print_insp_state.to_string() << std::endl;

    // muon
    const int pdg = interactor_state.pdg;

    // mass
    const scalar mass = interactor_state.mass;

    // new momentum
    const scalar newP = state._stepping._bound_params.charge() /
                        state._stepping._bound_params.qop();

    // new energy
    const scalar newE = std::sqrt(newP * newP + mass * mass);

    // Initial energy
    const scalar iniE = std::sqrt(iniP * iniP + mass * mass);

    // New qop variance
    const scalar new_var_qop =
        matrix_operator().element(state._stepping._bound_params.covariance(),
                                  e_bound_qoverp, e_bound_qoverp);

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // Same material used for default telescope detector
    material_slab<scalar> slab(mat, thickness);

    // Expected Bethe Stopping power for telescope geometry is estimated
    // as (number of planes * energy loss per plane assuming 1 GeV muon).
    // It is not perfectly precise as the track loses its energy during
    // propagation. However, since the energy loss << the track momentum,
    // the assumption is not very bad

    // -1 is required because the last surface is a portal
    const scalar dE =
        I.compute_energy_loss_bethe(is, slab, pdg, mass, q / iniP, q) *
        (positions.size() - 1);

    // Check if the new energy after propagation is enough close to the
    // expected value
    EXPECT_NEAR(newE, iniE - dE, 1e-5);

    const scalar sigma_qop = I.compute_energy_loss_landau_sigma_QOverP(
        is, slab, pdg, mass, q / iniP, q);

    const scalar dvar_qop = sigma_qop * sigma_qop * (positions.size() - 1);

    EXPECT_NEAR(new_var_qop, dvar_qop, 1e-10);

    // @todo: Validate the backward direction case as well?
}

scalar get_variance(const std::vector<scalar>& v) {
    scalar sum = std::accumulate(v.begin(), v.end(), 0.0);
    scalar mean = sum / v.size();
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   [mean](double x) { return x - mean; });
    scalar sq_sum =
        std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    scalar variance = sq_sum / v.size();
    return variance;
}

// Material interaction test with telescope Geometry
TEST(material_interaction, telescope_geometry_scattering_angle) {
    vecmem::host_memory_resource host_mr;

    using ln_stepper_t = line_stepper<transform3>;

    typename ln_stepper_t::free_track_parameters_type default_trk{
        {0, 0, 0}, 0, {1, 0, 0}, -1};

    // Build from given module positions
    std::vector<scalar> positions = {0., 1000. * unit_constants::cm};

    const auto mat = silicon_tml<scalar>();
    const scalar thickness = 500 * unit_constants::cm;
    // Use unbounded surfaces
    constexpr bool unbounded = true;

    const auto det = create_telescope_detector<unbounded>(
        host_mr, positions, ln_stepper_t(),
        typename ln_stepper_t::state{default_trk}, 2000. * unit_constants::mm,
        2000. * unit_constants::mm, mat, thickness);

    using navigator_t = navigator<decltype(det)>;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_default_policy;
    using stepper_t = line_stepper<transform3, constraints_t, policy_t>;
    using interactor_t = pointwise_material_interactor<transform3>;
    using simulator_t = scattering_simulator<interactor_t>;
    using actor_chain_t =
        actor_chain<dtuple, propagation::print_inspector, pathlimit_aborter,
                    parameter_transporter<transform3>, interactor_t,
                    simulator_t, resetter<transform3>>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator is built from the stepper and navigator
    propagator_t p({}, {});

    const scalar q = -1.;
    const scalar iniP = 10 * unit_constants::GeV;

    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = 0.;
    getter::element(bound_vector, e_bound_loc1, 0) = 0.;
    getter::element(bound_vector, e_bound_phi, 0) = 0;
    getter::element(bound_vector, e_bound_theta, 0) = M_PI_2;
    getter::element(bound_vector, e_bound_qoverp, 0) = q / iniP;
    getter::element(bound_vector, e_bound_time, 0) = 0.;
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    // bound track parameter
    const bound_track_parameters<transform3> bound_param(0, bound_vector,
                                                         bound_cov);

    int n_samples = 100000;
    std::vector<scalar> phi_vec;
    std::vector<scalar> theta_vec;

    scalar ref_phi_var(0);
    scalar ref_theta_var(0);

    for (int i = 0; i < n_samples; i++) {

        propagation::print_inspector::state print_insp_state{};
        pathlimit_aborter::state aborter_state{};
        parameter_transporter<transform3>::state bound_updater{};
        interactor_t::state interactor_state{};
        interactor_state.do_energy_loss = false;
        simulator_t::state simulator_state(interactor_state);
        resetter<transform3>::state resetter_state{};

        // Create actor states tuples
        actor_chain_t::state actor_states =
            std::tie(print_insp_state, aborter_state, bound_updater,
                     interactor_state, simulator_state, resetter_state);

        propagator_t::state state(bound_param, det, actor_states);

        state._stepping().set_overstep_tolerance(-1000. * unit_constants::um);

        // Propagate the entire detector
        ASSERT_TRUE(p.propagate(state))
            << print_insp_state.to_string() << std::endl;

        const auto& final_params = state._stepping._bound_params;

        if (i == 0) {
            const auto& covariance = final_params.covariance();
            ref_phi_var =
                matrix_operator().element(covariance, e_bound_phi, e_bound_phi);
            ref_theta_var = matrix_operator().element(covariance, e_bound_theta,
                                                      e_bound_theta);
        }

        phi_vec.push_back(final_params.phi());
        theta_vec.push_back(final_params.theta());
    }

    scalar phi_var = get_variance(phi_vec);
    scalar theta_var = get_variance(theta_vec);

    EXPECT_NEAR((phi_var - ref_phi_var) / ref_phi_var, 0, 0.05);
    EXPECT_NEAR((theta_var - ref_theta_var) / ref_theta_var, 0, 0.05);

    // To make sure that the varainces are zero
    EXPECT_TRUE(ref_phi_var > 1e-4 && ref_theta_var > 1e-4);
}
