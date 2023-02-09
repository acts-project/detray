/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/geometry/volume_graph.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/inspectors.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

// Vecmem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

// Type Definitions
template <typename mask_shape_t>
using tel_detector_t =
    detector<detector_registry::template telescope_detector<mask_shape_t>,
             covfie::field>;
/*using registry_type = detector_registry::default_detector;
using detector_t = detector<registry_type, covfie::field>;
using mask_container = typename detector_t::mask_container;
using material_container = typename detector_t::material_container;
using surface_type = typename detector_t::surface_container::value_type;
using mask_link_type = typename surface_type::mask_link;
using material_link_type = typename surface_type::material_link;*/
/*using mag_field_t = detector_t::bfield_type;

using matrix_operator = standard_matrix_operator<scalar>;
using transform3 = typename detector_t::transform3;
using vector2 = __plugin::vector2<scalar>;
using vector3 = __plugin::vector3<scalar>;

using navigator_t = navigator<detector_t>;
// Re-init the volume at every step, so that surface will be put into the
// navigation cache
using crk_stepper_t = rk_stepper<mag_field_t::view_t, transform3,
                                 constrained_step<>, always_init>;
using actor_chain_t =
    actor_chain<dtuple, target_aborter, parameter_transporter<transform3>,
                parameter_resetter<transform3>>;
using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain_t>;*/

using matrix_operator = standard_matrix_operator<scalar>;

namespace env {

constexpr scalar tol{1e-2f};
constexpr scalar rk_tolerance{1e-8f};

// Surface material
const material<scalar> mat = silicon<scalar>();
const scalar thickness{2.f * unit<scalar>::mm};

// memory resource
vecmem::host_memory_resource resource;

}  // namespace env

// This tests the path correction of cartesian coordinate
TEST(path_correction, cartesian2D) {

    // Create a detector
    /*detector_t det(env::resource, std::move(env::mag_field));

    // Mask and material ID
    constexpr auto mask_id = registry_type::mask_ids::e_rectangle2;
    constexpr auto material_id = registry_type::material_ids::e_slab;

    // Add a volume
    det.new_volume(
        volume_id::e_cylinder,
        {0.f, 0.f, 0.f, 0.f, -constant<scalar>::pi, constant<scalar>::pi});

    typename detector_t::surface_container surfaces(&env::resource);
    typename detector_t::transform_container transforms(env::resource);
    typename detector_t::mask_container masks(env::resource);
    typename detector_t::material_container materials(env::resource);

    // Add a surface
    mask_link_type mask_link{mask_id, 0u};
    material_link_type material_link{material_id, 0u};
    surfaces.emplace_back(0u, mask_link, material_link, 0u, dindex_invalid,
                          surface_id::e_sensitive);

    // Add a transform
    const vector3 t{0.f, 0.f, 0.f};
    const vector3 z{1.f, 0.f, 0.f};
    const vector3 x{0.f, 1.f, 0.f};
    transforms.emplace_back(env::ctx, t, z, x);

    // Add a mask
    const scalar hx{100.f * unit<scalar>::mm};
    const scalar hy{100.f * unit<scalar>::mm};
    masks.template emplace_back<mask_id>(empty_context{}, 0u, hx, hy);

    // Add a material
    const material<scalar> mat = silicon<scalar>();
    const scalar thickness{2.f * unit<scalar>::mm};
    materials.template emplace_back<material_id>(empty_context{}, mat,
                                                 thickness);

    typename detector_t::volume_type &vol = det.volume_by_index(0u);
    det.add_objects_per_volume(env::ctx, vol, surfaces, masks, transforms,
                               materials);*/
    using detector_t = tel_detector_t<unbounded<rectangle2D<>>>;
    using mag_field_t = detector_t::bfield_type;

    using transform3 = typename detector_t::transform3;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = __plugin::vector3<scalar>;

    using navigator_t = navigator<detector_t, navigation::print_inspector>;
    // Re-init the volume at every step, so that surface will be put into the
    // navigation cache
    using crk_stepper_t = rk_stepper<mag_field_t::view_t, transform3,
                                     constrained_step<>, always_init>;
    using actor_chain_t =
        actor_chain<dtuple, target_aborter, parameter_transporter<transform3>,
                    parameter_resetter<transform3>>;
    using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain_t>;

    // B field
    const vector3 B{0, 0, 1. * unit<scalar>::T};
    mag_field_t mag_field(
        mag_field_t::backend_t::configuration_t{B[0], B[1], B[2]});

    // Use unbounded rectangle surface
    mask<unbounded<rectangle2D<>>> rectangle{};

    // Build telescope detector with a single unbounded rectangle
    const auto det = create_telescope_detector(
        env::resource, std::move(mag_field), rectangle, 1u,
        500. * unit<scalar>::mm, env::mat, env::thickness);

    /// Prints linking information for every node when visited
    struct volume_printout {
        void operator()(const detector_t::volume_type &n) const {
            std::cout << "On volume: " << n.index() << std::endl;
        }
    };

    // Build the graph
    volume_graph<detector_t> graph(det);
    std::cout << graph.to_string() << std::endl;

    // Generate track starting point
    vector2 local{2.f, 3.f};
    vector3 mom{1.f * unit<scalar>::MeV, 0.f, 0.f};

    scalar time{0.f};
    scalar q{-1.f};

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0u) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0u) = local[1];
    getter::element(bound_vector, e_bound_phi, 0u) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0u) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0u) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0u) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.f;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.f;

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0u, bound_vector,
                                                          bound_cov);

    // Actors
    target_aborter::state targeter{0u};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, det.get_bfield(), det);

    crk_stepper_t::state &crk_state = propagation._stepping;

    // Path length per turn
    scalar S{2.f * getter::perp(mom) / getter::norm(env::B) *
             constant<scalar>::pi};
    // Constrain the step size to force frequent re-navigation for this
    // strongly bent track
    crk_state.template set_constraint<step::constraint::e_accuracy>(0.1f * S);

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = env::rk_tolerance;
    crk_state().set_overstep_tolerance(-1000000000000.f);

    // Run propagator
    // p.propagate(propagation, std::tie(targeter, bound_updater, rst));

    // Retrieve navigation information
    auto &debug_printer = propagation._navigation.inspector();

    ASSERT_TRUE(
        p.propagate(propagation, std::tie(targeter, bound_updater, rst)));
    std::cout << debug_printer.to_string();

    // Bound state after one turn propagation
    /*const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0u),
                    matrix_operator().element(bound_vec1, i, 0u), env::tol);
    }

    // covaraince
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), env::tol);
        }
    }

    /*
    // Print the matrix elements
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            printf("%f ", matrix_operator().element(bound_cov0, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            printf("%f ", matrix_operator().element(bound_cov1, i, j));
        }
        printf("\n");
    }
    printf("\n");
    */
}

/*TEST(path_correction, polar) {

    // Create a detector
    detector_t det(env::resource);

    // Mask and material ID
    constexpr auto mask_id = registry_type::mask_ids::e_ring2;
    constexpr auto material_id = registry_type::material_ids::e_slab;

    // Add a volume
    det.new_volume(
        volume_id::e_cylinder,
        {0.f, 0.f, 0.f, 0.f, -constant<scalar>::pi, constant<scalar>::pi});

    typename detector_t::surface_container surfaces(&env::resource);
    typename detector_t::transform_container transforms(env::resource);
    typename detector_t::mask_container masks(env::resource);
    typename detector_t::material_container materials(env::resource);

    // Add a surface
    mask_link_type mask_link{mask_id, 0u};
    material_link_type material_link{material_id, 0u};
    surfaces.emplace_back(0u, mask_link, material_link, 0, dindex_invalid,
                          surface_id::e_sensitive);

    // Add a transform
    const vector3 t{0.f, 0.f, 0.f};
    const vector3 z{1.f, 0.f, 0.f};
    const vector3 x{0.f, 1.f, 0.f};
    transforms.emplace_back(env::ctx, t, z, x);

    // Add a mask
    const scalar r_low{0.f * unit<scalar>::mm};
    const scalar r_high{100.f * unit<scalar>::mm};
    masks.template emplace_back<mask_id>(empty_context{}, 0u, r_low, r_high);

    // Add a material
    const material<scalar> mat = silicon<scalar>();
    const scalar thickness{2.f * unit<scalar>::mm};
    materials.template emplace_back<material_id>(empty_context{}, mat,
                                                 thickness);

    typename detector_t::volume_type &vol = det.volume_by_index(0);
    det.add_objects_per_volume(env::ctx, vol, surfaces, masks, transforms,
                               materials);

    // Generate track starting point
    vector2 local{2.f, constant<scalar>::pi / 6.f};
    vector3 mom{1.f * unit<scalar>::MeV, 0.f, 0.f};

    scalar time{0.f};
    scalar q{-1.f};

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0u) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0u) = local[1];
    getter::element(bound_vector, e_bound_phi, 0u) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0u) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0u) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0u) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.f;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.f;

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0u, bound_vector,
                                                          bound_cov);

    // Actors
    target_aborter::state targeter{0u};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, env::mag_field, det);

    crk_stepper_t::state &crk_state = propagation._stepping;

    // RK stepper and its state
    EXPECT_NEAR(crk_state().pos()[0], 0.f, env::tol);
    EXPECT_NEAR(crk_state().pos()[1],
                2.f * std::cos(constant<scalar>::pi / 6.f), env::tol);
    EXPECT_NEAR(crk_state().pos()[2],
                2.f * std::sin(constant<scalar>::pi / 6.f), env::tol);

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = env::rk_tolerance;

    // Path length per turn
    scalar S{2.f * getter::perp(mom) / getter::norm(env::B) *
             constant<scalar>::pi};
    // Constrain the step size to force frequent re-navigation for this
    // strongly bent track
    crk_state.template set_constraint<step::constraint::e_accuracy>(0.1f * S);

    // Run propagator
    p.propagate(propagation, std::tie(targeter, bound_updater, rst));

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0u),
                    matrix_operator().element(bound_vec1, i, 0u), env::tol);
    }

    // covaraince
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), env::tol);
        }
    }
}

TEST(path_correction, cylindrical) {

    // Create a detector
    detector_t det(env::resource);

    // Mask and material ID
    constexpr auto mask_id = registry_type::mask_ids::e_cylinder2;
    constexpr auto material_id = registry_type::material_ids::e_slab;

    // Add a volume
    det.new_volume(
        volume_id::e_cylinder,
        {0.f, 0.f, 0.f, 0.f, -constant<scalar>::pi, constant<scalar>::pi});

    typename detector_t::surface_container surfaces(&env::resource);
    typename detector_t::transform_container transforms(env::resource);
    typename detector_t::mask_container masks(env::resource);
    typename detector_t::material_container materials(env::resource);

    // Add a surface
    mask_link_type mask_link{mask_id, 0u};
    material_link_type material_link{material_id, 0u};
    surfaces.emplace_back(0u, mask_link, material_link, 0u, dindex_invalid,
                          surface_id::e_sensitive);

    // Add a transform
    const vector3 t{-50.f * unit<scalar>::mm, 0.f, 0.f};
    const vector3 z{0.f, 0.f, 1.f};
    const vector3 x{1.f, 0.f, 0.f};
    transforms.emplace_back(env::ctx, t, z, x);

    // Add a mask
    const scalar r{50.f * unit<scalar>::mm};
    const scalar half_length_1{1000.f * unit<scalar>::mm};
    const scalar half_length_2{1000.f * unit<scalar>::mm};
    masks.template emplace_back<mask_id>(empty_context{}, 0u, r, half_length_1,
                                         half_length_2);

    // Add a material
    const material<scalar> mat = silicon<scalar>();
    const scalar thickness{2.f * unit<scalar>::mm};
    materials.template emplace_back<material_id>(empty_context{}, mat,
                                                 thickness);

    typename detector_t::volume_type &vol = det.volume_by_index(0);
    det.add_objects_per_volume(env::ctx, vol, surfaces, masks, transforms,
                               materials);

    // Generate track starting point
    vector2 local{0.f, 0.f};
    vector3 mom{1.f * unit<scalar>::MeV, 0.f, 0.f};

    scalar time{0.f};
    scalar q{-1.f};

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0u) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0u) = local[1];
    getter::element(bound_vector, e_bound_phi, 0u) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0u) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0u) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0u) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.f;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.f;

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    // Actors
    target_aborter::state targeter{0u};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, env::mag_field, det);

    crk_stepper_t::state &crk_state = propagation._stepping;

    // RK stepper and its state
    EXPECT_NEAR(crk_state().pos()[0], 0.f, env::tol);  // x
    EXPECT_NEAR(crk_state().pos()[1], 0.f, env::tol);  // y
    EXPECT_NEAR(crk_state().pos()[2], 0.f, env::tol);  // z

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = env::rk_tolerance;

    // Path length per turn
    scalar S{2.f * getter::perp(mom) / getter::norm(env::B) *
             constant<scalar>::pi};
    // Constrain the step size to force frequent re-navigation for this
    // strongly bent track
    crk_state.template set_constraint<step::constraint::e_accuracy>(0.1f * S);

    // Run propagator
    p.propagate(propagation, std::tie(targeter, bound_updater, rst));

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0u),
                    matrix_operator().element(bound_vec1, i, 0u), env::tol);
    }

    // covaraince
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), env::tol);
        }
    }
}*/