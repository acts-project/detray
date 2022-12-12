/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

// Vecmem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;

struct surface_targeter : actor {

    struct state {
        scalar _path;
        dindex _target_surface_index = dindex_invalid;
    };

    /// Enforces thepath limit on a stepper state
    ///
    /// @param abrt_state contains the path limit
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &actor_state,
                                       propagator_state_t &prop_state) const {

        prop_state._heartbeat = true;
        auto &navigation = prop_state._navigation;
        auto &stepping = prop_state._stepping;
        navigation.set_full_trust();

        scalar residual = actor_state._path - stepping.path_length();
        stepping.set_constraint(residual);

        typename propagator_state_t::navigator_state_type::intersection_t is;
        is.index = actor_state._target_surface_index;
        is.path = residual;
        auto &candidates = navigation.candidates();
        candidates.clear();
        candidates.push_back(is);
        navigation.set_next(candidates.begin());
        navigation.set_unknown();

        if (residual < std::abs(1e-6)) {
            prop_state._heartbeat = false;
            candidates.push_back({});
            navigation.next()++;
            navigation.set_on_module();
        }
    }
};

// Type Definitions
using vector2 = __plugin::vector2<scalar>;
using vector3 = __plugin::vector3<scalar>;
using matrix_operator = standard_matrix_operator<scalar>;
using registry_type = detector_registry::default_detector;
using detector_type = detector<registry_type, covfie::field>;
using mask_container = typename detector_type::mask_container;
using material_container = typename detector_type::material_container;
using transform3 = typename detector_type::transform3;
using surface_type = typename detector_type::surface_container::value_type;
using mask_link_type = typename surface_type::mask_link;
using material_link_type = typename surface_type::material_link;
using mag_field_t = detector_type::bfield_type;
using navigator_t = navigator<detector_type>;
using crk_stepper_t =
    rk_stepper<mag_field_t::view_t, transform3, constrained_step<>>;
using actor_chain_t =
    actor_chain<dtuple, surface_targeter, parameter_transporter<transform3>,
                parameter_resetter<transform3>>;
using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain_t>;

namespace env {

constexpr scalar epsilon = 1e-2;
constexpr scalar rk_tolerance = 1e-8;

// B field
const vector3 B{0, 0, 1. * unit<scalar>::T};
mag_field_t mag_field(mag_field_t::backend_t::configuration_t{B[0], B[1],
                                                              B[2]});

// memory resource
vecmem::host_memory_resource resource;

// Context
const typename detector_type::geometry_context ctx{};

}  // namespace env

// This tests the path correction of cartesian coordinate
TEST(path_correction, cartesian) {

    // Create a detector
    detector_type det(env::resource, std::move(env::mag_field));

    // Mask and material ID
    constexpr auto mask_id = registry_type::mask_ids::e_rectangle2;
    constexpr auto material_id = registry_type::material_ids::e_slab;

    // Add a volume
    det.new_volume(volume_id::e_cylinder, {0., 0., 0., 0., -M_PI, M_PI});

    typename detector_type::surface_container surfaces(&env::resource);
    typename detector_type::transform_container transforms(env::resource);
    typename detector_type::mask_container masks(env::resource);
    typename detector_type::material_container materials(env::resource);

    // Add a surface
    const auto trf_index = transforms.size(env::ctx);
    mask_link_type mask_link{mask_id, masks.template size<mask_id>()};
    material_link_type material_link{material_id,
                                     materials.template size<material_id>()};
    surfaces.emplace_back(trf_index, mask_link, material_link, 0,
                          dindex_invalid, surface_id::e_sensitive);

    // Add a transform
    const vector3 t{0, 0, 0};
    const vector3 z{1, 0, 0};
    const vector3 x{0, 1, 0};
    transforms.emplace_back(env::ctx, t, z, x);

    // Add a mask
    const scalar hx = 100. * unit<scalar>::mm;
    const scalar hy = 100. * unit<scalar>::mm;
    masks.template emplace_back<mask_id>(empty_context{}, 0UL, hx, hy);

    // Add a material
    const material<scalar> mat = silicon<scalar>();
    const scalar thickness = 2 * unit<scalar>::mm;
    materials.template emplace_back<material_id>(empty_context{}, mat,
                                                 thickness);

    typename detector_type::volume_type &vol = det.volume_by_index(0);
    det.add_objects_per_volume(env::ctx, vol, surfaces, masks, transforms,
                               materials);

    // Generate track starting point
    vector2 local{2, 3};
    vector3 mom{1 * unit<scalar>::MeV, 0., 0.};

    scalar time = 0.;
    scalar q = -1.;

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0) = local[1];
    getter::element(bound_vector, e_bound_phi, 0) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.;

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    // Path length per turn
    scalar S = 2. * getter::perp(mom) / getter::norm(env::B) * M_PI;

    // Actors
    surface_targeter::state targeter{S, 0};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, env::mag_field, det);

    crk_stepper_t::state &crk_state = propagation._stepping;

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = env::rk_tolerance;

    // Run propagator
    p.propagate(propagation, std::tie(targeter, bound_updater, rst));

    EXPECT_NEAR(crk_state.path_length(), targeter._path, env::epsilon);

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (std::size_t i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0),
                    matrix_operator().element(bound_vec1, i, 0), env::epsilon);
    }

    // covaraince
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j),
                        env::epsilon);
        }
    }

    /*
    // Print the matrix elements
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            printf("%f ", matrix_operator().element(bound_cov0, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            printf("%f ", matrix_operator().element(bound_cov1, i, j));
        }
        printf("\n");
    }
    printf("\n");
    */
}

TEST(path_correction, polar) {

    // Create a detector
    detector_type det(env::resource);

    // Mask and material ID
    constexpr auto mask_id = registry_type::mask_ids::e_ring2;
    constexpr auto material_id = registry_type::material_ids::e_slab;

    // Add a volume
    det.new_volume(volume_id::e_cylinder, {0., 0., 0., 0., -M_PI, M_PI});

    typename detector_type::surface_container surfaces(&env::resource);
    typename detector_type::transform_container transforms(env::resource);
    typename detector_type::mask_container masks(env::resource);
    typename detector_type::material_container materials(env::resource);

    // Add a surface
    const auto trf_index = transforms.size(env::ctx);
    mask_link_type mask_link{mask_id, masks.template size<mask_id>()};
    material_link_type material_link{material_id,
                                     materials.template size<material_id>()};
    surfaces.emplace_back(trf_index, mask_link, material_link, 0,
                          dindex_invalid, surface_id::e_sensitive);

    // Add a transform
    const vector3 t{0, 0, 0};
    const vector3 z{1, 0, 0};
    const vector3 x{0, 1, 0};
    transforms.emplace_back(env::ctx, t, z, x);

    // Add a mask
    const scalar r_low = 0. * unit<scalar>::mm;
    const scalar r_high = 100. * unit<scalar>::mm;
    masks.template emplace_back<mask_id>(empty_context{}, 0UL, r_low, r_high);

    // Add a material
    const material<scalar> mat = silicon<scalar>();
    const scalar thickness = 2 * unit<scalar>::mm;
    materials.template emplace_back<material_id>(empty_context{}, mat,
                                                 thickness);

    typename detector_type::volume_type &vol = det.volume_by_index(0);
    det.add_objects_per_volume(env::ctx, vol, surfaces, masks, transforms,
                               materials);

    // Generate track starting point
    vector2 local{2, M_PI / 6.};
    vector3 mom{1 * unit<scalar>::MeV, 0., 0.};

    scalar time = 0.;
    scalar q = -1.;

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0) = local[1];
    getter::element(bound_vector, e_bound_phi, 0) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.;

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    // Path length per turn
    scalar S = 2. * getter::perp(mom) / getter::norm(env::B) * M_PI;

    // Actors
    surface_targeter::state targeter{S, 0};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, env::mag_field, det);

    crk_stepper_t::state &crk_state = propagation._stepping;

    // RK stepper and its state
    EXPECT_FLOAT_EQ(crk_state().pos()[0], 0);
    EXPECT_FLOAT_EQ(crk_state().pos()[1], 2 * std::cos(M_PI / 6.));
    EXPECT_FLOAT_EQ(crk_state().pos()[2], 2 * std::sin(M_PI / 6.));

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = env::rk_tolerance;

    // Run propagator
    p.propagate(propagation, std::tie(targeter, bound_updater, rst));

    EXPECT_NEAR(crk_state.path_length(), targeter._path, env::epsilon);

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (std::size_t i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0),
                    matrix_operator().element(bound_vec1, i, 0), env::epsilon);
    }

    // covaraince
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j),
                        env::epsilon);
        }
    }
}

TEST(path_correction, cylindrical) {

    // Create a detector
    detector_type det(env::resource);

    // Mask and material ID
    constexpr auto mask_id = registry_type::mask_ids::e_cylinder2;
    constexpr auto material_id = registry_type::material_ids::e_slab;

    // Add a volume
    det.new_volume(volume_id::e_cylinder, {0., 0., 0., 0., -M_PI, M_PI});

    typename detector_type::surface_container surfaces(&env::resource);
    typename detector_type::transform_container transforms(env::resource);
    typename detector_type::mask_container masks(env::resource);
    typename detector_type::material_container materials(env::resource);

    // Add a surface
    const auto trf_index = transforms.size(env::ctx);
    mask_link_type mask_link{mask_id, masks.template size<mask_id>()};
    material_link_type material_link{material_id,
                                     materials.template size<material_id>()};
    surfaces.emplace_back(trf_index, mask_link, material_link, 0,
                          dindex_invalid, surface_id::e_sensitive);

    // Add a transform
    const vector3 t{-50 * unit<scalar>::mm, 0, 0};
    const vector3 z{0, 0, 1};
    const vector3 x{1, 0, 0};
    transforms.emplace_back(env::ctx, t, z, x);

    // Add a mask
    const scalar r = 50 * unit<scalar>::mm;
    const scalar half_length_1 = 1000. * unit<scalar>::mm;
    const scalar half_length_2 = 1000. * unit<scalar>::mm;
    masks.template emplace_back<mask_id>(empty_context{}, 0UL, r, half_length_1,
                                         half_length_2);

    // Add a material
    const material<scalar> mat = silicon<scalar>();
    const scalar thickness = 2 * unit<scalar>::mm;
    materials.template emplace_back<material_id>(empty_context{}, mat,
                                                 thickness);

    typename detector_type::volume_type &vol = det.volume_by_index(0);
    det.add_objects_per_volume(env::ctx, vol, surfaces, masks, transforms,
                               materials);

    // Generate track starting point
    vector2 local{0., 0.};
    vector3 mom{1 * unit<scalar>::MeV, 0., 0.};
    scalar time = 0.;
    scalar q = -1.;

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0) = local[1];
    getter::element(bound_vector, e_bound_phi, 0) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.;

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    // Path length per turn
    scalar S = 2. * getter::perp(mom) / getter::norm(env::B) * M_PI;

    // Actors
    surface_targeter::state targeter{S, 0};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, env::mag_field, det);

    crk_stepper_t::state &crk_state = propagation._stepping;

    // RK stepper and its state
    EXPECT_FLOAT_EQ(crk_state().pos()[0], 0);
    EXPECT_FLOAT_EQ(crk_state().pos()[1], 0);  // y
    EXPECT_FLOAT_EQ(crk_state().pos()[2], 0);  // z

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = env::rk_tolerance;

    // Run propagator
    p.propagate(propagation, std::tie(targeter, bound_updater, rst));

    EXPECT_NEAR(crk_state.path_length(), targeter._path, env::epsilon);

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (std::size_t i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0),
                    matrix_operator().element(bound_vec1, i, 0), env::epsilon);
    }

    // covaraince
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j),
                        env::epsilon);
        }
    }
}