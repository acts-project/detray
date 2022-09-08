/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/inspectors.hpp"

// VecMem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using point2 = __plugin::point2<scalar>;
using point3 = __plugin::point3<scalar>;
using vector3 = __plugin::vector3<scalar>;
using transform3 = __plugin::transform3<scalar>;

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
        I.compute_energy_loss_bethe(is, slab, pdg, m, qOverP) /
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
    const scalar dE = I.compute_energy_loss_landau(is, slab, pdg, m, qOverP) /
                      (unit_constants::MeV);

    // Check if difference is within 5% error
    EXPECT_TRUE(std::abs(std::get<2>(GetParam()) - dE) / dE < 0.05);

    // Landau Energy loss Fluctuation
    const scalar fwhm =
        I.compute_energy_loss_landau_fwhm(is, slab, pdg, m, qOverP) /
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
TEST(material_interaction, telescope_geometry) {

    vecmem::host_memory_resource host_mr;

    using ln_stepper_t = line_stepper<transform3>;

    typename ln_stepper_t::free_track_parameters_type default_trk{
        {0, 0, 0}, 0, {0, 0, 1}, -1};

    // Build from given module positions
    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};

    const auto det =
        create_telescope_detector(host_mr, positions, ln_stepper_t(),
                                  typename ln_stepper_t::state{default_trk});

    using navigator_t = navigator<decltype(det)>;
    using constraints_t = constrained_step<>;
    using policy_t = stepper_default_policy;
    using stepper_t = line_stepper<transform3, constraints_t, policy_t>;
    using interactor_t = pointwise_material_interactor<interaction<scalar>>;
    using actor_chain_t = actor_chain<dtuple, propagation::print_inspector,
                                      pathlimit_aborter, interactor_t>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator is built from the stepper and navigator
    propagator_t p({}, {});

    const scalar q = -1.;
    const point3 pos = {0., 0., 0.};
    const vector3 mom = {0., 0., 1.};
    free_track_parameters<transform3> track(pos, 0, mom, q);

    propagation::print_inspector::state print_insp_state{};
    pathlimit_aborter::state aborter_state{};
    interactor_t::state interactor_state{};

    // Create actor states tuples
    actor_chain_t::state actor_states =
        std::tie(print_insp_state, aborter_state, interactor_state);

    propagator_t::state state(track, det, actor_states);

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

    // Initial momentum
    const scalar iniP = getter::norm(mom);

    // Initial energy
    const scalar iniE = std::sqrt(iniP * iniP + mass * mass);

    // Interaction object
    interaction<scalar> I;

    // intersection with a zero incidence angle
    line_plane_intersection is;

    // Same material used for default telescope detector
    material_slab<scalar> slab(silicon_tml<scalar>(), 80 * unit_constants::um);

    // Expected Bethe Stopping power for telescope geometry is estimated
    // as (number of planes * energy loss per plane assuming 1 GeV muon).
    // It is not perfectly precise as the track loses its energy during
    // propagation. However, since the energy loss << the track momentum,
    // the assumption is not very bad
    const scalar dE =
        I.compute_energy_loss_bethe(is, slab, pdg, mass, q / iniP) *
        (positions.size() -
         1);  // -1 is required because the last surface is a portal

    // Check if the new energy after propagation is enough close to the
    // expected value
    EXPECT_NEAR(newE, iniE - dE, 1e-5);

    // @todo: Validate the backward direction case as well?
}
