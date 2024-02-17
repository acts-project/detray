/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/navigation/intersection/helix_intersector.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/utils/inspectors.hpp"
#include "detray/utils/statistics.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// Boost
#include <boost/program_options.hpp>

// System include(s).
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace po = boost::program_options;
using namespace detray;

// Type declarations
using transform3_type = __plugin::transform3<scalar>;
using vector3 = __plugin::vector3<scalar>;
using bound_vector_type = bound_track_parameters<transform3_type>::vector_type;
using bound_covariance_type =
    bound_track_parameters<transform3_type>::covariance_type;
using matrix_operator = typename transform3_type::matrix_actor;
using size_type = typename transform3_type::size_type;
template <size_type ROWS, size_type COLS>
using matrix_type = typename matrix_operator::template matrix_type<ROWS, COLS>;

namespace {

// B-field Setup
// 1.996 T from averaged Bz strength of ODD field in the region of sqrt(x^2 +
// y^2) < 500 mm and abs(z) < 500 mm
const vector3 B_z{0.f, 0.f, 1.996f * unit<scalar>::T};

// Initial delta for numerical differentiaion
const std::array<scalar, 5u> h_sizes_rect{1e-1f, 1e-1f, 1e-2f, 1e-3f, 1e-2f};
const std::array<scalar, 5u> h_sizes_wire{1e-1f, 1e-1f, 1e-2f, 1e-3f, 1e-2f};

// Ridders' algorithm setup
constexpr const unsigned int Nt = 10u;
const std::array<scalar, 5u> safe{2.0f, 2.0f, 2.0f, 2.0f, 2.0f};
const std::array<scalar, 5u> con{1.2f, 1.2f, 1.2f, 1.2f, 1.2f};
constexpr const scalar big = std::numeric_limits<scalar>::max();

std::random_device rd;
// For detector generation
std::mt19937_64 mt1(rd());
// For smearing initial parameter
std::mt19937_64 mt2(rd());

// Detector length generator
constexpr const scalar min_detector_length = 50.f * unit<scalar>::mm;
constexpr const scalar max_detector_length = 500.f * unit<scalar>::mm;
std::uniform_real_distribution<scalar> rand_length(min_detector_length,
                                                   max_detector_length);

// Euler angles for the surface rotation
std::uniform_real_distribution<scalar> rand_alpha(0.f,
                                                  2.f * constant<scalar>::pi);
std::uniform_real_distribution<scalar> rand_cosbeta(0.5f, 1.f);
std::uniform_int_distribution<int> rand_bool(0, 1);
std::uniform_real_distribution<scalar> rand_gamma(0.f,
                                                  2.f * constant<scalar>::pi);

// constexpr const scalar shift = 5.f * unit<scalar>::mm;
//  Shift is disabled at the moment
constexpr const scalar shift = 0.f * unit<scalar>::mm;
std::uniform_real_distribution<scalar> rand_shift(-shift, shift);

// surface types
using rect_type = rectangle2D;
using wire_type = line<true>;

}  // namespace

void wrap_angles(const bound_vector_type& ref_vector,
                 bound_vector_type& target_vector) {

    const scalar rphi = getter::element(ref_vector, e_bound_phi, 0u);
    const scalar tphi = getter::element(target_vector, e_bound_phi, 0u);
    scalar new_tphi = tphi;

    if (rphi >= constant<scalar>::pi_2) {
        if (tphi < 0.f) {
            new_tphi += 2.f * constant<scalar>::pi;
        }
    } else if (rphi < -constant<scalar>::pi_2) {
        if (tphi >= 0.f) {
            new_tphi -= 2.f * constant<scalar>::pi;
        }
    }

    getter::element(target_vector, e_bound_phi, 0u) = new_tphi;
}

scalar get_relative_difference(scalar ref_val, scalar num_val) {
    scalar rel_diff{0.f};

    // If the evaulated jacovian or numerical diffentiation is too small set the
    // relative difference to zero
    if (ref_val == 0.f || num_val == 0.f) {
        rel_diff = 0.f;
    } else {
        rel_diff = std::abs(ref_val - num_val) / std::abs(num_val);
    }

    return rel_diff;
}

// Get random initial covariance
bound_covariance_type get_random_initial_covariance(const scalar ini_qop) {

    // Initial covariance matrix for smearing
    bound_covariance_type ini_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();

    // Correlation factor in the range of [-10%, 10%]
    scalar min_corr = -0.1f;
    scalar max_corr = 0.1f;
    // Random correction factor
    std::uniform_real_distribution<scalar> rand_corr(min_corr, max_corr);

    std::normal_distribution<scalar> rand_l0(0.f * unit<scalar>::um,
                                             50.f * unit<scalar>::um);
    std::normal_distribution<scalar> rand_l1(0.f * unit<scalar>::um,
                                             50.f * unit<scalar>::um);
    std::normal_distribution<scalar> rand_phi(0.f * unit<scalar>::mrad,
                                              1.f * unit<scalar>::mrad);
    std::normal_distribution<scalar> rand_theta(0.f * unit<scalar>::mrad,
                                                1.f * unit<scalar>::mrad);
    std::normal_distribution<scalar> rand_qop(0.f, 0.01f * ini_qop);
    std::normal_distribution<scalar> rand_time(0.f * unit<scalar>::ns,
                                               1.f * unit<scalar>::ns);

    /*
    // Typical stddev range taken from the figures of ATL-PHYS-PUB-2021-024 and
    // ATLAS-TDR-030
    std::normal_distribution<scalar> rand_l0(5.f * unit<scalar>::um,
                                             200.f * unit<scalar>::um);
    std::normal_distribution<scalar> rand_l1(10.f * unit<scalar>::um,
                                             4000.f * unit<scalar>::um);
    std::normal_distribution<scalar> rand_phi(0.05f * unit<scalar>::mrad,
                                              5.0f * unit<scalar>::mrad);
    std::normal_distribution<scalar> rand_theta(0.01f * unit<scalar>::mrad,
                                                2.0f * unit<scalar>::mrad);
    std::normal_distribution<scalar> rand_qop(0.01f * ini_qop, 0.1f * ini_qop);
    std::normal_distribution<scalar> rand_time(0.f * unit<scalar>::ns,
                                               1.f * unit<scalar>::ns);
    */

    std::array<scalar, 6u> stddevs;
    stddevs[0] = rand_l0(mt2);
    stddevs[1] = rand_l1(mt2);
    stddevs[2] = rand_phi(mt2);
    stddevs[3] = rand_theta(mt2);
    stddevs[4] = rand_qop(mt2);
    stddevs[5] = rand_time(mt2);

    for (unsigned int i = 0u; i < 6u; i++) {
        for (unsigned int j = 0u; j < 6u; j++) {
            if (i == j) {
                getter::element(ini_cov, i, i) = stddevs[i] * stddevs[i];
            } else if (i > j) {
                getter::element(ini_cov, i, j) =
                    stddevs[i] * stddevs[j] * rand_corr(mt2);
                getter::element(ini_cov, j, i) = getter::element(ini_cov, i, j);
            }
        }
    }

    return ini_cov;
}

// Input covariance should be the diagonal matrix
bound_vector_type get_smeared_bound_vector(const bound_covariance_type& ini_cov,
                                           const bound_vector_type& ini_vec) {

    // Do the Cholesky Decomposition
    const bound_covariance_type L =
        matrix_helper<matrix_operator>().cholesky_decompose(ini_cov);

    // Vecor with random elements from a normal distribution
    bound_vector_type k = matrix_operator().template zero<e_bound_size, 1u>();
    std::normal_distribution<scalar> normal_dist(0.f, 1.f);
    for (unsigned int i = 0u; i < 5u; i++) {
        // Smear the value
        getter::element(k, i, 0) = normal_dist(mt2);
    }

    const bound_vector_type new_vec = ini_vec + L * k;

    return new_vec;
}

template <typename detector_t, typename detector_t::metadata::mask_ids mask_id>
std::pair<euler_rotation<transform3_type>, std::array<scalar, 3u>> tilt_surface(
    detector_t& det, const unsigned int sf_id, const vector3& helix_dir) {

    const auto& sf = det.surface(sf_id);
    const auto& trf_link = sf.transform();
    auto& trf = det.transform_store()[trf_link];

    euler_rotation<transform3_type> euler;
    euler.alpha = rand_alpha(mt1);

    if (sf_id == 1u) {
        const int is_counterclockwise = rand_bool(mt1);
        auto beta = math::acos(rand_cosbeta(mt1));
        if (is_counterclockwise == 0) {
            beta = -beta;
        }
        euler.beta = beta;
        euler.gamma = rand_gamma(mt1);
    }

    // Helix direction
    euler.z = helix_dir;

    if constexpr (mask_id == detector_t::masks::id::e_rectangle2) {
        // ubasis = trf.x() for bound frame
        euler.x = trf.x();
    } else if (mask_id == detector_t::masks::id::e_cell_wire) {
        // ubasis = vxt/|vxt| where v is trf.z() and t is helix_dir
        EXPECT_NEAR(vector::dot(trf.z(), helix_dir), 0.f, 1e-6f);
        euler.x = vector::normalize(vector::cross(trf.z(), helix_dir));
    }

    EXPECT_NEAR(vector::dot(euler.x, euler.z), 0.f, 1e-6f);

    auto [local_x, local_z] = euler();

    if constexpr (mask_id == detector_t::masks::id::e_cell_wire) {
        // local_z should be the line direction
        local_z = vector::cross(local_z, local_x);
    }

    scalar x_shift{0.f};
    scalar y_shift{0.f};
    scalar z_shift{0.f};

    if (sf_id == 1u) {
        x_shift = rand_shift(mt1);
        y_shift = rand_shift(mt1);
        z_shift = rand_shift(mt1);
    }
    // Translation vector
    vector3 translation = trf.translation();

    vector3 displacement({x_shift, y_shift, z_shift});
    translation = trf.translation() + displacement;
    transform3_type new_trf(translation, local_z, local_x, true);
    det.transform_store()[trf_link] = new_trf;

    return {euler, {x_shift, y_shift, z_shift}};
}

template <typename transform3_t>
struct bound_getter : actor {

    // Transformation matching this struct
    using transform3_type = transform3_t;
    // scalar_type
    using scalar_type = typename transform3_type::scalar_type;
    using bound_track_parameters_type = bound_track_parameters<transform3_t>;
    using free_track_parameters_type = free_track_parameters<transform3_t>;

    struct state {

        scalar m_min_path_length;
        scalar m_path_length;
        bound_track_parameters_type m_param_departure;
        bound_track_parameters_type m_param_destination;
        typename bound_track_parameters_type::covariance_type m_jacobi;
        scalar m_avg_step_size{0.f};
        std::size_t step_count{0u};
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& actor_state,
                                       propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        actor_state.step_count++;

        const scalar N = static_cast<scalar>(actor_state.step_count);

        actor_state.m_avg_step_size =
            ((N - 1.f) * actor_state.m_avg_step_size + stepping._step_size) / N;

        if (navigation.is_on_module() && navigation.barcode().index() == 0u) {

            actor_state.m_param_departure = stepping._bound_params;
        }
        // Get the bound track parameters and jacobian at the destination
        // surface
        else if (navigation.is_on_module() &&
                 navigation.barcode().index() == 1u) {

            actor_state.m_path_length = stepping._path_length;
            actor_state.m_param_destination = stepping._bound_params;
            actor_state.m_jacobi = stepping._full_jacobian;

            // Stop navigation if the destination surface found
            propagation._heartbeat &= navigation.exit();
        }

        if (stepping._path_length > actor_state.m_min_path_length) {
            propagation._navigation.set_no_trust();
        }

        return;
    }
};

/// Numerically integrate the jacobian
template <typename propagator_t, typename field_t>
bound_getter<transform3_type>::state evaluate_bound_param(
    const scalar detector_length,
    const bound_track_parameters<transform3_type>& initial_param,
    const typename propagator_t::detector_type& det, const field_t& field,
    const scalar overstep_tolerance, const scalar on_surface_tolerance,
    const scalar rk_tolerance, const scalar constraint_step,
    bool use_field_gradient, bool do_scatter, bool do_covariance_transport,
    bool do_inspect = false) {

    // Propagator is built from the stepper and navigator
    propagation::config<scalar> cfg{};
    cfg.navigation.overstep_tolerance = overstep_tolerance;
    cfg.navigation.on_surface_tolerance = on_surface_tolerance;
    cfg.stepping.rk_error_tol = rk_tolerance;
    cfg.stepping.use_eloss_gradient = true;
    cfg.stepping.use_field_gradient = use_field_gradient;
    cfg.stepping.do_scatter = do_scatter;
    cfg.stepping.do_covariance_transport = do_covariance_transport;
    propagator_t p(cfg);

    // Actor states
    parameter_transporter<transform3_type>::state transporter_state{};
    bound_getter<transform3_type>::state bound_getter_state{};
    bound_getter_state.m_min_path_length = detector_length * 0.75f;
    parameter_resetter<transform3_type>::state resetter_state{};
    auto actor_states =
        std::tie(transporter_state, bound_getter_state, resetter_state);

    // Init propagator states for the reference track
    typename propagator_t::state state(initial_param, field, det);

    // Run the propagation for the reference track
    state.do_debug = do_inspect;
    state._stepping
        .template set_constraint<detray::step::constraint::e_accuracy>(
            constraint_step);

    p.propagate(state, actor_states);
    if (do_inspect) {
        std::cout << state.debug_stream.str() << std::endl;
    }

    return bound_getter_state;
}

template <typename propagator_t, typename field_t>
typename bound_track_parameters<transform3_type>::vector_type
get_displaced_bound_vector(
    const bound_track_parameters<transform3_type>& ref_param,
    const typename propagator_t::detector_type& det,
    const scalar detector_length, const field_t& field,
    const scalar overstep_tolerance, const scalar on_surface_tolerance,
    const scalar rk_tolerance, const scalar constraint_step,
    const unsigned int target_index, const scalar displacement) {

    propagation::config<scalar> cfg{};
    cfg.navigation.overstep_tolerance = overstep_tolerance;
    cfg.navigation.on_surface_tolerance = on_surface_tolerance;
    cfg.stepping.rk_error_tol = rk_tolerance;
    // we don't have to scatering or covariance transport
    cfg.stepping.do_scatter = false;
    cfg.stepping.do_covariance_transport = false;

    // Propagator is built from the stepper and navigator
    propagator_t p(cfg);

    bound_track_parameters<transform3_type> dparam = ref_param;
    auto dvec = dparam.vector();
    getter::element(dvec, target_index, 0u) += displacement;

    dparam.set_vector(dvec);

    typename propagator_t::state dstate(dparam, field, det);

    // Actor states
    parameter_transporter<transform3_type>::state transporter_state{};
    parameter_resetter<transform3_type>::state resetter_state{};
    bound_getter<transform3_type>::state bound_getter_state{};
    bound_getter_state.m_min_path_length = detector_length * 0.75f;

    auto actor_states =
        std::tie(transporter_state, bound_getter_state, resetter_state);
    dstate._stepping
        .template set_constraint<detray::step::constraint::e_accuracy>(
            constraint_step);

    p.propagate(dstate, actor_states);

    auto new_vec = bound_getter_state.m_param_destination.vector();

    // phi needs to be wrapped w.r.t. phi of the reference vector
    wrap_angles(ref_param.vector(), new_vec);

    return new_vec;
}

/// Numerically evaluate the jacobian
template <typename propagator_t, typename field_t>
bound_track_parameters<transform3_type>::covariance_type directly_differentiate(
    const bound_track_parameters<transform3_type>& ref_param,
    const typename propagator_t::detector_type& det,
    const scalar detector_length, const field_t& field,
    const scalar overstep_tolerance, const scalar on_surface_tolerance,
    const scalar rk_tolerance, const scalar constraint_step,
    const std::array<scalar, 5u> hs,
    [[maybe_unused]] const scalar avg_step_size,
    std::array<bool, 25>& convergence) {

    // Return Jacobian
    bound_track_parameters<transform3_type>::covariance_type
        differentiated_jacobian;

    for (unsigned int i = 0u; i < 5u; i++) {

        scalar delta = hs[i];

        if (i == e_bound_theta) {
            const scalar rtheta = ref_param.theta();
            if (rtheta < constant<scalar>::pi_2) {
                delta = math::min(delta, rtheta);
            } else if (rtheta >= constant<scalar>::pi_2) {
                delta = math::min(delta, constant<scalar>::pi - rtheta);
            }
        }

        // Containers for Ridders algorithm
        // Page 231 of [Numerical Recipes] 3rd edition
        std::array<std::array<std::array<scalar, Nt>, Nt>, 5u> Arr;

        scalar fac;
        std::array<scalar, 5u> errt;
        std::array<scalar, 5u> err{big, big, big, big, big};
        std::array<bool, 5u> complete{false, false, false, false, false};

        const auto vec1 = get_displaced_bound_vector<propagator_t, field_t>(
            ref_param, det, detector_length, field, overstep_tolerance,
            on_surface_tolerance, rk_tolerance, constraint_step, i,
            1.f * delta);
        const auto vec2 = get_displaced_bound_vector<propagator_t, field_t>(
            ref_param, det, detector_length, field, overstep_tolerance,
            on_surface_tolerance, rk_tolerance, constraint_step, i,
            -1.f * delta);

        for (unsigned int j = 0; j < 5u; j++) {

            const scalar v1 = getter::element(vec1, j, 0u);
            const scalar v2 = getter::element(vec2, j, 0u);

            Arr[j][0][0] = (v1 - v2) / (2.f * delta);
        }

        for (unsigned int p = 1u; p < Nt; p++) {
            delta /= con[i];

            const auto nvec1 =
                get_displaced_bound_vector<propagator_t, field_t>(
                    ref_param, det, detector_length, field, overstep_tolerance,
                    on_surface_tolerance, rk_tolerance, constraint_step, i,
                    1.f * delta);
            const auto nvec2 =
                get_displaced_bound_vector<propagator_t, field_t>(
                    ref_param, det, detector_length, field, overstep_tolerance,
                    on_surface_tolerance, rk_tolerance, constraint_step, i,
                    -1.f * delta);

            for (unsigned int j = 0; j < 5u; j++) {

                const scalar v1 = getter::element(nvec1, j, 0u);
                const scalar v2 = getter::element(nvec2, j, 0u);

                Arr[j][0][p] = (v1 - v2) / (2.f * delta);
            }

            const scalar con2 = con[i] * con[i];
            fac = con2;

            for (unsigned int q = 1; q <= p; q++) {

                for (unsigned int j = 0; j < 5u; j++) {
                    Arr[j][q][p] =
                        (Arr[j][q - 1][p] * fac - Arr[j][q - 1][p - 1]) /
                        (fac - 1.0f);
                    fac = con2 * fac;

                    errt[j] = math::max(
                        math::abs(Arr[j][q][p] - Arr[j][q - 1][p]),
                        math::abs(Arr[j][q][p] - Arr[j][q - 1][p - 1]));

                    if (errt[j] <= err[j] && complete[j] == false) {
                        err[j] = errt[j];
                        getter::element(differentiated_jacobian, j, i) =
                            Arr[j][q][p];
                    }
                }
            }

            for (unsigned int j = 0; j < 5u; j++) {
                if (math::abs(Arr[j][p][p] - Arr[j][p - 1][p - 1]) >=
                    safe[i] * err[j]) {
                    complete[j] = true;
                }
            }

            if (std::count(complete.begin(), complete.end(), false) == 0u) {
                break;
            }
        }

        // Row-major
        for (std::size_t j = 0u; j < 5u; j++) {
            convergence[i + j * 5] = complete[j];
        }

        /*
        if (std::count(complete.begin(), complete.end(), false) > 0u) {
            std::cout << "The jacobian is not converged!" << std::endl;
        } else {
            std::cout << "The jacobian is converged!" << std::endl;
        }
        */
    }

    return differentiated_jacobian;
}

template <typename detector_t, typename detector_t::metadata::mask_ids mask_id>
bound_track_parameters<transform3_type> get_initial_parameter(
    detector_t& det, const free_track_parameters<transform3_type>& vertex,
    const vector3& field, const scalar helix_tolerance) {

    // Helix from the vertex
    detail::helix<transform3_type> hlx(vertex, &field);

    const auto& departure_sf = det.surface(0u);
    const auto& trf_link = departure_sf.transform();
    const auto& departure_trf = det.transform_store()[trf_link];
    const auto& mask_link = departure_sf.mask();
    const auto& departure_mask =
        det.mask_store().template get<mask_id>().at(mask_link.index());

    using mask_t =
        typename detector_t::mask_container::template get_type<mask_id>;
    helix_intersector<typename mask_t::shape, transform3_type> hlx_is{};
    hlx_is.convergence_tolerance = helix_tolerance;
    auto sfi = hlx_is(hlx, departure_sf, departure_mask, departure_trf, 0.f);
    EXPECT_EQ(sfi.status, intersection::status::e_inside)
        << " Initial surface not found" << std::endl
        << " log10(Helix tolerance): " << math::log10(helix_tolerance)
        << " Phi: " << getter::phi(vertex.dir())
        << " Theta: " << getter::theta(vertex.dir())
        << " Mom [GeV/c]: " << vertex.p() << std::endl
        << sfi;

    const auto path_length = sfi.path;
    // As we don't rotate or shift the initial surface anymore, the path_length
    // should be 0
    EXPECT_FLOAT_EQ(static_cast<float>(path_length), 0.f);

    const auto pos = hlx(path_length);
    const auto dir = hlx.dir(path_length);

    const free_track_parameters<transform3_type> free_par(pos, 0, dir,
                                                          hlx._qop);

    const auto bound_vec =
        surface{det, departure_sf}.free_to_bound_vector({}, free_par.vector());

    bound_track_parameters<transform3_type> ret;
    ret.set_surface_link(geometry::barcode{0u});
    ret.set_vector(bound_vec);

    return ret;
}

template <typename propagator_t, typename field_t>
void evaluate_jacobian_difference(
    const unsigned int trk_count, typename propagator_t::detector_type& det,
    const scalar detector_length,
    const bound_track_parameters<transform3_type>& track, const field_t& field,
    const material<scalar> volume_mat, const scalar overstep_tolerance,
    const scalar on_surface_tolerance, const scalar rk_tolerance,
    const scalar constraint_step, const std::array<scalar, 5u>& hs,
    std::ofstream& file, scalar& ref_rel_diff, bool use_field_gradient = true,
    bool do_inspect = false) {

    const auto phi0 = track.phi();
    const auto theta0 = track.theta();
    (void)phi0;
    (void)theta0;

    det.volumes()[0u].set_material(volume_mat);

    const bool do_scatter = true;
    const bool do_covariance_transport = true;

    // No scattering for jacobian comparsion study
    // But perform covariance transport
    auto bound_getter = evaluate_bound_param<propagator_t, field_t>(
        detector_length, track, det, field, overstep_tolerance,
        on_surface_tolerance, rk_tolerance, constraint_step, use_field_gradient,
        !do_scatter, do_covariance_transport, do_inspect);

    const auto reference_param = bound_getter.m_param_departure;
    const auto final_param = bound_getter.m_param_destination;

    // Sanity check
    ASSERT_EQ(reference_param.surface_link().index(), 0u)
        << " Initial surface not found " << std::endl
        << " log10(RK tolerance): " << math::log10(rk_tolerance)
        << " Path length [mm]: " << bound_getter.m_path_length
        << " Average step size [mm]: " << bound_getter.m_avg_step_size
        << " Phi: " << reference_param.phi()
        << " Theta: " << reference_param.theta()
        << " Mom [GeV/c]: " << reference_param.p();
    ASSERT_EQ(final_param.surface_link().index(), 1u)
        << " Final surface not found " << std::endl
        << " log10(RK tolerance): " << math::log10(rk_tolerance)
        << " Path length [mm]: " << bound_getter.m_path_length
        << " Average step size [mm]: " << bound_getter.m_avg_step_size
        << " Phi: " << reference_param.phi()
        << " Theta: " << reference_param.theta()
        << " Mom [GeV/c]: " << reference_param.p();
    ASSERT_GE(bound_getter.m_path_length, 0.f);
    ASSERT_LE(bound_getter.m_path_length, max_detector_length + 100.f);

    const auto reference_jacobian = bound_getter.m_jacobi;

    file << trk_count << ",";

    file << reference_param.bound_local()[0] << ","
         << reference_param.bound_local()[1] << "," << reference_param.phi()
         << "," << reference_param.theta() << "," << reference_param.qop()
         << ",";

    file << final_param.bound_local()[0] << "," << final_param.bound_local()[1]
         << "," << final_param.phi() << "," << final_param.theta() << ","
         << final_param.qop() << ",";

    std::array<bool, 25u> convergence;

    auto differentiated_jacobian =
        directly_differentiate<propagator_t, field_t>(
            reference_param, det, detector_length, field, overstep_tolerance,
            on_surface_tolerance, rk_tolerance, constraint_step, hs,
            bound_getter.m_avg_step_size, convergence);

    bool total_convergence =
        (std::count(convergence.begin(), convergence.end(), false) == 0);

    // Convergence
    file << total_convergence << ",";
    for (unsigned int i = 0; i < 25u; i++) {
        file << convergence[i] << ",";
    }

    // Reference track
    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {
            file << getter::element(reference_jacobian, i, j) << ",";
        }
    }

    // Numerical evaluation
    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {
            file << getter::element(differentiated_jacobian, i, j) << ",";
        }
    }

    // Difference between evaluation and direct jacobian
    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {

            const scalar ref_val = getter::element(reference_jacobian, i, j);
            const scalar num_val =
                getter::element(differentiated_jacobian, i, j);
            const scalar rel_diff = get_relative_difference(ref_val, num_val);

            file << rel_diff << ",";

            // We return dqopdqop for test
            if (i == 4 && j == 4) {
                ref_rel_diff = rel_diff;
            }
        }
    }

    // Path length
    file << bound_getter.m_path_length << ",";

    // Average step size
    file << bound_getter.m_avg_step_size << ",";

    // Log10(RK tolerance)
    file << math::log10(rk_tolerance) << ",";

    // Log10(on surface tolerance)
    file << math::log10(on_surface_tolerance) << ",";

    // Overstep tolerance
    file << overstep_tolerance;

    file << std::endl;
}

template <typename propagator_t, typename field_t>
void evaluate_covariance_transport(
    const unsigned int trk_count, typename propagator_t::detector_type& det,
    const scalar detector_length,
    const bound_track_parameters<transform3_type>& track, const field_t& field,
    const material<scalar> volume_mat, const scalar overstep_tolerance,
    const scalar on_surface_tolerance, const scalar rk_tolerance,
    const scalar constraint_step, std::ofstream& file,
    bool use_field_gradient = true) {

    det.volumes()[0u].set_material(volume_mat);

    // Copy track
    auto track_copy = track;

    // Make initial covariance
    const bound_covariance_type ini_cov =
        get_random_initial_covariance(track_copy.qop());

    track_copy.set_covariance(ini_cov);

    const bool do_scatter = true;
    const bool do_covariance_transport = true;

    // Reference track (Without Scattering but with covariance transport)
    auto bound_getter = evaluate_bound_param<propagator_t, field_t>(
        detector_length, track_copy, det, field, overstep_tolerance,
        on_surface_tolerance, rk_tolerance, constraint_step, use_field_gradient,
        !do_scatter, do_covariance_transport);

    const auto reference_param = bound_getter.m_param_departure;
    const auto ini_vec = reference_param.vector();
    const auto final_param = bound_getter.m_param_destination;
    const auto fin_vec = final_param.vector();
    const auto fin_cov = final_param.covariance();

    // Sanity check
    ASSERT_EQ(reference_param.surface_link().index(), 0u)
        << " Initial surface not found " << std::endl
        << " log10(RK tolerance): " << math::log10(rk_tolerance)
        << " Path length [mm]: " << bound_getter.m_path_length
        << " Average step size [mm]: " << bound_getter.m_avg_step_size
        << " Phi: " << reference_param.phi()
        << " Theta: " << reference_param.theta()
        << " Mom [GeV/c]: " << reference_param.p();
    ASSERT_EQ(final_param.surface_link().index(), 1u)
        << " Final surface not found " << std::endl
        << " log10(RK tolerance): " << math::log10(rk_tolerance)
        << " Path length [mm]: " << bound_getter.m_path_length
        << " Average step size [mm]: " << bound_getter.m_avg_step_size
        << " Phi: " << reference_param.phi()
        << " Theta: " << reference_param.theta()
        << " Mom [GeV/c]: " << reference_param.p();
    ASSERT_GE(bound_getter.m_path_length, min_detector_length - 10.f * shift);
    ASSERT_LE(bound_getter.m_path_length, max_detector_length + 10.f * shift);

    // Get smeared initial bound vector
    const bound_vector_type smeared_ini_vec =
        get_smeared_bound_vector(ini_cov, reference_param.vector());

    // Make smeared bound track parameter
    auto smeared_track = track_copy;
    smeared_track.set_vector(smeared_ini_vec);

    // Smeared track (With scattering and no covariance transport)
    auto smeared_bound_getter = evaluate_bound_param<propagator_t, field_t>(
        detector_length, smeared_track, det, field, overstep_tolerance,
        on_surface_tolerance, rk_tolerance, constraint_step, use_field_gradient,
        do_scatter, !do_covariance_transport);

    // Get smeared final bound vector
    bound_vector_type smeared_fin_vec =
        smeared_bound_getter.m_param_destination.vector();

    // phi needs to be wrapped w.r.t. phi of the reference vector
    wrap_angles(fin_vec, smeared_fin_vec);

    // Get pull values
    std::array<scalar, 5u> pulls;

    const bound_vector_type diff = smeared_fin_vec - fin_vec;
    for (unsigned int i = 0u; i < 5u; i++) {
        pulls[i] = getter::element(diff, i, 0u) /
                   math::sqrt(getter::element(fin_cov, i, i));
    }

    // Get Chi2
    const matrix_type<1u, 1u> chi2 = matrix_operator().transpose(diff) *
                                     matrix_operator().inverse(fin_cov) * diff;
    const scalar chi2_val = getter::element(chi2, 0u, 0u);

    file << trk_count << ",";

    // File writing
    file << getter::element(ini_vec, e_bound_loc0, 0u) << ","
         << getter::element(ini_vec, e_bound_loc1, 0u) << ","
         << getter::element(ini_vec, e_bound_phi, 0u) << ","
         << getter::element(ini_vec, e_bound_theta, 0u) << ","
         << getter::element(ini_vec, e_bound_qoverp, 0u) << ",";

    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {
            file << getter::element(ini_cov, i, j) << ",";
        }
    }

    file << getter::element(fin_vec, e_bound_loc0, 0u) << ","
         << getter::element(fin_vec, e_bound_loc1, 0u) << ","
         << getter::element(fin_vec, e_bound_phi, 0u) << ","
         << getter::element(fin_vec, e_bound_theta, 0u) << ","
         << getter::element(fin_vec, e_bound_qoverp, 0u) << ",";

    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {
            file << getter::element(fin_cov, i, j) << ",";
        }
    }

    file << getter::element(smeared_ini_vec, e_bound_loc0, 0u) << ","
         << getter::element(smeared_ini_vec, e_bound_loc1, 0u) << ","
         << getter::element(smeared_ini_vec, e_bound_phi, 0u) << ","
         << getter::element(smeared_ini_vec, e_bound_theta, 0u) << ","
         << getter::element(smeared_ini_vec, e_bound_qoverp, 0u) << ",";

    file << getter::element(smeared_fin_vec, e_bound_loc0, 0u) << ","
         << getter::element(smeared_fin_vec, e_bound_loc1, 0u) << ","
         << getter::element(smeared_fin_vec, e_bound_phi, 0u) << ","
         << getter::element(smeared_fin_vec, e_bound_theta, 0u) << ","
         << getter::element(smeared_fin_vec, e_bound_qoverp, 0u) << ",";

    file << pulls[0] << "," << pulls[1] << "," << pulls[2] << "," << pulls[3]
         << "," << pulls[4] << ",";

    file << chi2_val << ",";

    // Path length
    file << bound_getter.m_path_length << ",";

    // Average step size
    file << bound_getter.m_avg_step_size << ",";

    // Log10(RK tolerance)
    file << math::log10(rk_tolerance) << ",";

    // Log10(on surface tolerance)
    file << math::log10(on_surface_tolerance) << ",";

    // Overstep tolerance
    file << overstep_tolerance;

    file << std::endl;
}

template <typename detector_t, typename detector_t::metadata::mask_ids mask_id>
typename bound_track_parameters<transform3_type>::vector_type
get_displaced_bound_vector_helix(
    const bound_track_parameters<transform3_type>& track, const vector3& field,
    unsigned int target_index, scalar displacement, const detector_t& det,
    const scalar helix_tolerance) {

    const auto& departure_sf = det.surface(0u);

    const auto& destination_sf = det.surface(1u);
    const auto& trf_link = destination_sf.transform();
    const auto& destination_trf = det.transform_store()[trf_link];
    const auto& mask_link = destination_sf.mask();
    const auto& destination_mask =
        det.mask_store().template get<mask_id>().at(mask_link.index());

    auto dvec = track.vector();
    getter::element(dvec, target_index, 0u) += displacement;
    const auto free_vec =
        surface{det, departure_sf}.bound_to_free_vector({}, dvec);
    detail::helix<transform3_type> hlx(free_vec, &field);

    using mask_t =
        typename detector_t::mask_container::template get_type<mask_id>;
    helix_intersector<typename mask_t::shape, transform3_type> hlx_is{};
    hlx_is.convergence_tolerance = helix_tolerance;
    auto sfi =
        hlx_is(hlx, destination_sf, destination_mask, destination_trf, 0.f);
    const auto path_length = sfi.path;
    const auto pos = hlx(path_length);
    const auto dir = hlx.dir(path_length);

    const free_track_parameters<transform3_type> new_free_par(pos, 0, dir,
                                                              hlx._qop);
    auto new_bound_vec = surface{det, destination_sf}.free_to_bound_vector(
        {}, new_free_par.vector());

    // phi needs to be wrapped w.r.t. phi of the reference vector
    wrap_angles(dvec, new_bound_vec);

    return new_bound_vec;
}

template <typename detector_t, typename detector_t::metadata::mask_ids mask_id>
void evaluate_jacobian_difference_helix(
    const unsigned int trk_count, detector_t& det,
    const bound_track_parameters<transform3_type>& track, const vector3& field,
    const std::array<scalar, 5u> hs, std::ofstream& file,
    const scalar helix_tolerance) {

    const auto phi0 = track.phi();
    const auto theta0 = track.theta();
    (void)phi0;
    (void)theta0;

    // Get bound to free Jacobi
    const auto& departure_sf = det.surface(0u);
    const auto bound_to_free_jacobi =
        surface{det, departure_sf}.bound_to_free_jacobian({}, track.vector());

    // Get fre vector
    const auto free_vec =
        surface{det, departure_sf}.bound_to_free_vector({}, track.vector());
    // Helix from the departure surface
    detail::helix<transform3_type> hlx(free_vec, &field);

    const auto& destination_sf = det.surface(1u);
    const auto& trf_link = destination_sf.transform();
    const auto& destination_trf = det.transform_store()[trf_link];
    const auto& mask_link = destination_sf.mask();
    const auto& destination_mask =
        det.mask_store().template get<mask_id>().at(mask_link.index());

    using mask_t =
        typename detector_t::mask_container::template get_type<mask_id>;
    helix_intersector<typename mask_t::shape, transform3_type> hlx_is{};
    hlx_is.convergence_tolerance = helix_tolerance;

    auto sfi =
        hlx_is(hlx, destination_sf, destination_mask, destination_trf, 0.f);

    EXPECT_EQ(sfi.status, intersection::status::e_inside)
        << " Final surface not found" << std::endl
        << " log10(Helix tolerance): " << math::log10(helix_tolerance)
        << " Phi: " << track.phi() << " Theta: " << track.theta()
        << " Mom [GeV/c]: " << track.p();

    const auto path_length = sfi.path;

    // Get transport Jacobi
    const auto transport_jacobi = hlx.jacobian(path_length);

    const auto pos = hlx(path_length);
    const auto dir = hlx.dir(path_length);
    const auto qop = hlx._qop;

    // Get correction term
    const auto correction_term =
        matrix_operator().template identity<e_free_size, e_free_size>() +
        surface{det, destination_sf}.path_correction(
            {}, pos, dir, qop * vector::cross(dir, field), 0.f);

    const free_track_parameters<transform3_type> free_par(pos, 0.f, dir, qop);

    // Get free to bound Jacobi
    const auto free_to_bound_jacobi =
        surface{det, destination_sf}.free_to_bound_jacobian({},
                                                            free_par.vector());

    // Get full Jacobi
    const auto reference_jacobian = free_to_bound_jacobi * correction_term *
                                    transport_jacobi * bound_to_free_jacobi;

    // Get bound vector
    const auto bound_vec = surface{det, destination_sf}.free_to_bound_vector(
        {}, free_par.vector());

    /******************************
     *  Numerical differentiation
     * ****************************/

    bound_track_parameters<transform3_type>::covariance_type
        differentiated_jacobian;

    std::array<bool, 25u> convergence;

    for (unsigned int i = 0u; i < 5u; i++) {

        scalar delta = hs[i];

        if (i == e_bound_theta) {
            const scalar rtheta = track.theta();
            if (rtheta < constant<scalar>::pi_2) {
                delta = math::min(delta, rtheta);
            } else if (rtheta >= constant<scalar>::pi_2) {
                delta = math::min(delta, constant<scalar>::pi - rtheta);
            }
        }

        // Containers for Ridders algorithm
        // Page 231 of [Numerical Recipes] 3rd edition
        std::array<std::array<std::array<scalar, Nt>, Nt>, 5u> Arr;

        scalar fac;
        std::array<scalar, 5u> errt;
        std::array<scalar, 5u> err{big, big, big, big, big};
        std::array<bool, 5u> complete{false, false, false, false, false};

        const auto vec1 = get_displaced_bound_vector_helix<detector_t, mask_id>(
            track, field, i, 1.f * delta, det, helix_tolerance);
        const auto vec2 = get_displaced_bound_vector_helix<detector_t, mask_id>(
            track, field, i, -1.f * delta, det, helix_tolerance);

        for (unsigned int j = 0; j < 5u; j++) {

            const scalar v1 = getter::element(vec1, j, 0u);
            const scalar v2 = getter::element(vec2, j, 0u);

            Arr[j][0][0] = (v1 - v2) / (2.f * delta);
        }

        for (unsigned int p = 1u; p < Nt; p++) {
            delta /= con[i];

            const auto nvec1 =
                get_displaced_bound_vector_helix<detector_t, mask_id>(
                    track, field, i, 1.f * delta, det, helix_tolerance);
            const auto nvec2 =
                get_displaced_bound_vector_helix<detector_t, mask_id>(
                    track, field, i, -1.f * delta, det, helix_tolerance);

            for (unsigned int j = 0; j < 5u; j++) {

                const scalar v1 = getter::element(nvec1, j, 0u);
                const scalar v2 = getter::element(nvec2, j, 0u);

                Arr[j][0][p] = (v1 - v2) / (2.f * delta);
            }

            const scalar con2 = con[i] * con[i];
            fac = con2;

            for (unsigned int q = 1; q <= p; q++) {

                for (unsigned int j = 0; j < 5u; j++) {
                    Arr[j][q][p] =
                        (Arr[j][q - 1][p] * fac - Arr[j][q - 1][p - 1]) /
                        (fac - 1.0f);
                    fac = con2 * fac;

                    errt[j] = math::max(
                        math::abs(Arr[j][q][p] - Arr[j][q - 1][p]),
                        math::abs(Arr[j][q][p] - Arr[j][q - 1][p - 1]));

                    if (errt[j] <= err[j] && complete[j] == false) {
                        err[j] = errt[j];
                        getter::element(differentiated_jacobian, j, i) =
                            Arr[j][q][p];
                    }
                }
            }

            for (unsigned int j = 0; j < 5u; j++) {

                if (math::abs(Arr[j][p][p] - Arr[j][p - 1][p - 1]) >=
                    safe[i] * err[j]) {
                    complete[j] = true;
                }
            }

            if (std::count(complete.begin(), complete.end(), false) == 0u) {
                break;
            }
        }

        for (std::size_t j = 0u; j < 5u; j++) {
            convergence[i + j * 5] = complete[j];
        }
    }

    bool total_convergence =
        (std::count(convergence.begin(), convergence.end(), false) == 0);

    file << trk_count << ",";

    file << track.bound_local()[0] << "," << track.bound_local()[1] << ","
         << track.phi() << "," << track.theta() << "," << track.qop() << ",";

    file << getter::element(bound_vec, e_bound_loc0, 0u) << ","
         << getter::element(bound_vec, e_bound_loc1, 0u) << ","
         << getter::element(bound_vec, e_bound_phi, 0u) << ","
         << getter::element(bound_vec, e_bound_theta, 0u) << ","
         << getter::element(bound_vec, e_bound_qoverp, 0u) << ",";

    // Convergence
    file << total_convergence << ",";
    for (unsigned int i = 0; i < 25u; i++) {
        file << convergence[i] << ",";
    }

    // Reference track
    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {
            file << getter::element(reference_jacobian, i, j) << ",";
        }
    }

    // Numerical evaluation
    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {
            file << getter::element(differentiated_jacobian, i, j) << ",";
        }
    }

    // Difference between evaluation and direct jacobian
    for (unsigned int i = 0; i < 5u; i++) {
        for (unsigned int j = 0; j < 5u; j++) {

            const scalar ref_val = getter::element(reference_jacobian, i, j);
            const scalar num_val =
                getter::element(differentiated_jacobian, i, j);
            const scalar rel_diff = get_relative_difference(ref_val, num_val);
            file << rel_diff << ",";
        }
    }

    // Path length
    file << path_length << ",";

    // Average step size (Doesn't exist for helix intersection)
    file << 0 << ",";

    // Log10(RK tolerance) (Doesn't exist for helix intersection)
    file << 0 << ",";

    // Log10(helix intersection tolerance)
    file << math::log10(helix_tolerance) << ",";

    // Overstep tolerance (Doesn't exist for helix intersection)
    file << 0;

    file << std::endl;
}

void setup_csv_header_jacobian(std::ofstream& file) {

    file << std::fixed << std::showpoint;
    file << std::setprecision(32);

    // Track ID
    file << "track_ID,";

    // Initial Parameter at the departure surface
    file << "l0_I,l1_I,phi_I,theta_I,qop_I,";

    // Final Parameter at the destination surface
    file << "l0_F,l1_F,phi_F,theta_F,qop_F,";

    // Convergence
    file << "total_convergence,";
    file << "dl0dl0_C,dl0dl1_C,dl0dphi_C,dl0dtheta_C,dl0dqop_C,";
    file << "dl1dl0_C,dl1dl1_C,dl1dphi_C,dl1dtheta_C,dl1dqop_C,";
    file << "dphidl0_C,dphidl1_C,dphidphi_C,dphidtheta_C,dphidqop_C,";
    file << "dthetadl0_C,dthetadl1_C,dthetadphi_C,dthetadtheta_C,"
            "dthetadqop_C,";
    file << "dqopdl0_C,dqopdl1_C,dqopdphi_C,dqopdtheta_C,dqopdqop_C,";

    // Evaluation
    file << "dl0dl0_E,dl0dl1_E,dl0dphi_E,dl0dtheta_E,dl0dqop_E,";
    file << "dl1dl0_E,dl1dl1_E,dl1dphi_E,dl1dtheta_E,dl1dqop_E,";
    file << "dphidl0_E,dphidl1_E,dphidphi_E,dphidtheta_E,dphidqop_E,";
    file << "dthetadl0_E,dthetadl1_E,dthetadphi_E,dthetadtheta_E,"
            "dthetadqop_E,";
    file << "dqopdl0_E,dqopdl1_E,dqopdphi_E,dqopdtheta_E,dqopdqop_E,";

    // Numerical Differentiation
    file << "dl0dl0_D,dl0dl1_D,dl0dphi_D,dl0dtheta_D,dl0dqop_D,";
    file << "dl1dl0_D,dl1dl1_D,dl1dphi_D,dl1dtheta_D,dl1dqop_D,";
    file << "dphidl0_D,dphidl1_D,dphidphi_D,dphidtheta_D,dphidqop_D,";
    file << "dthetadl0_D,dthetadl1_D,dthetadphi_D,dthetadtheta_D,"
            "dthetadqop_D,";
    file << "dqopdl0_D,dqopdl1_D,dqopdphi_D,dqopdtheta_D,dqopdqop_D,";

    // Relative Difference between Evaluation and Numerical Differentiation
    file << "dl0dl0_R,dl0dl1_R,dl0dphi_R,dl0dtheta_R,dl0dqop_R,";
    file << "dl1dl0_R,dl1dl1_R,dl1dphi_R,dl1dtheta_R,dl1dqop_R,";
    file << "dphidl0_R,dphidl1_R,dphidphi_R,dphidtheta_R,dphidqop_R,";
    file << "dthetadl0_R,dthetadl1_R,dthetadphi_R,dthetadtheta_R,"
            "dthetadqop_R,";
    file << "dqopdl0_R,dqopdl1_R,dqopdphi_R,dqopdtheta_R,dqopdqop_R,";

    // Path length [mm]
    file << "path_length,";

    // Average step size [mm]
    file << "average_step_size,";

    // RK Tolerances [mm]
    file << "log10_rk_tolerance,";

    // Intersection tolerance [mm]
    file << "log10_intersection_tolerance,";

    // Overstep tolerance [mm]
    file << "overstep_tolerance";

    file << std::endl;
}

void setup_csv_header_covariance(std::ofstream& file) {

    file << std::fixed << std::showpoint;
    file << std::setprecision(32);

    // Track ID
    file << "track_ID,";

    // Initial parameters (vector + covariance) of reference track
    file << "l0_I,l1_I,phi_I,theta_I,qop_I,";
    file << "l0l0_I,l0l1_I,l0phi_I,l0theta_I,l0qop_I,";
    file << "l1l0_I,l1l1_I,l1phi_I,l1theta_I,l1qop_I,";
    file << "phil0_I,phil1_I,phiphi_I,phitheta_I,phiqop_I,";
    file << "thetal0_I,thetal1_I,thetaphi_I,thetatheta_I,thetaqop_I,";
    file << "qopl0_I,qopl1_I,qopphi_I,qoptheta_I,qopqop_I,";

    // Final parameters (vector + covariance) of reference track
    file << "l0_F,l1_F,phi_F,theta_F,qop_F,";
    file << "l0l0_F,l0l1_F,l0phi_F,l0theta_F,l0qop_F,";
    file << "l1l0_F,l1l1_F,l1phi_F,l1theta_F,l1qop_F,";
    file << "phil0_F,phil1_F,phiphi_F,phitheta_F,phiqop_F,";
    file << "thetal0_F,thetal1_F,thetaphi_F,thetatheta_F,thetaqop_F,";
    file << "qopl0_F,qopl1_F,qopphi_F,qoptheta_F,qopqop_F,";

    // Initial parameter (vector only) of smeared track
    file << "l0_IS,l1_IS,phi_IS,theta_IS,qop_IS,";

    // Final parameter (vector only) of smeared track
    file << "l0_FS,l1_FS,phi_FS,theta_FS,qop_FS,";

    // Pull values
    file << "pull_l0,pull_l1,pull_phi,pull_theta,pull_qop,";

    // Chi2
    file << "chi2,";

    // Path length [mm]
    file << "path_length,";

    // Average step size [mm]
    file << "average_step_size,";

    // Tolerances [mm]
    file << "log10_rk_tolerance,";

    // Intersection tolerance [mm]
    file << "log10_intersection_tolerance,";

    // Overstep tolerance [mm]
    file << "overstep_tolerance";

    file << std::endl;
}

int main(int argc, char** argv) {

    // Options parsing
    po::options_description desc("\ndetray jacobian validation options");
    desc.add_options()("help", "produce help message");
    desc.add_options()("output-directory",
                       po::value<std::string>()->default_value(""),
                       "Output directory");
    desc.add_options()("n-tracks",
                       po::value<std::size_t>()->default_value(100u),
                       "Number of tracks for generator");
    desc.add_options()("n-skips", po::value<std::size_t>()->default_value(0u),
                       "Number of skipped indices");
    desc.add_options()("skip-rect", po::value<bool>()->default_value(false),
                       "Skip rectangular telescope");
    desc.add_options()("skip-wire", po::value<bool>()->default_value(false),
                       "Skip wire telescope");
    desc.add_options()("log10-rk-tolerance-mm",
                       po::value<scalar>()->default_value(-4.f),
                       "Set log10(rk_tolerance_in_mm)");
    desc.add_options()("log10-helix-tolerance-mm",
                       po::value<scalar>()->default_value(-3.f),
                       "Set log10(helix_tolerance_in_mm)");
    desc.add_options()("overstep-tolerance-mm",
                       po::value<scalar>()->default_value(-1000.f),
                       "Set the overstep tolerance in mm unit");
    desc.add_options()("log10-on-surface-tolerance-mm",
                       po::value<scalar>()->default_value(-3.f),
                       "Set log10(on_surface_tolerance_in_mm)");
    desc.add_options()("rk-tolerance-iterate-mode",
                       po::value<bool>()->default_value(true),
                       "Iterate over the rk tolerances");
    desc.add_options()("log10-min-rk-tolerance-mm",
                       po::value<scalar>()->default_value(-6.f),
                       "Set log10(min_rk_tolerance_in_mm)");
    desc.add_options()("log10-max-rk-tolerance-mm",
                       po::value<scalar>()->default_value(2.f),
                       "Set log10(max_rk_tolerance_in_mm)");
    desc.add_options()("mc-seed", po::value<std::size_t>()->default_value(0u),
                       "Monte-Carlo seed");
    desc.add_options()("verbose-level", po::value<int>()->default_value(1),
                       "Verbose level");

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, desc,
                                 po::command_line_style::unix_style ^
                                     po::command_line_style::allow_short),
              vm);
    po::notify(vm);

    // Help message
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return EXIT_FAILURE;
    }

    const std::string output_directory =
        vm["output-directory"].as<std::string>();
    std::size_t n_tracks = vm["n-tracks"].as<std::size_t>();
    std::size_t n_skips = vm["n-skips"].as<std::size_t>();
    const bool skip_rect = vm["skip-rect"].as<bool>();
    const bool skip_wire = vm["skip-wire"].as<bool>();
    const scalar rk_power = vm["log10-rk-tolerance-mm"].as<scalar>();
    const scalar rk_tol = std::pow(10.f, rk_power) * unit<scalar>::mm;
    const scalar helix_power = vm["log10-helix-tolerance-mm"].as<scalar>();
    const scalar helix_tol = std::pow(10.f, helix_power) * unit<scalar>::mm;
    const scalar on_surface_power =
        vm["log10-on-surface-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
    const scalar on_surface_tol = math::pow(10.f, on_surface_power);
    const bool rk_tolerance_iterate_mode =
        vm["rk-tolerance-iterate-mode"].as<bool>();
    const scalar log10_min_rk_tolerance =
        vm["log10-min-rk-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
    const scalar log10_max_rk_tolerance =
        vm["log10-max-rk-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
    const std::size_t mc_seed = vm["mc-seed"].as<std::size_t>();
    const int verbose_lvl = vm["verbose-level"].as<int>();

    std::vector<scalar> log10_tols;
    scalar r = log10_min_rk_tolerance;
    while (r <= log10_max_rk_tolerance + 1e-3f) {
        log10_tols.push_back(r);
        r = r + 2.f;
    }

    // Set seed for random generator
    mt1.seed(mc_seed);

    // Volume material
    // const material<scalar> volume_mat = detray::silicon<scalar>();
    // const material<scalar> volume_mat = detray::iron_with_ded<scalar>();
    const material<scalar> volume_mat =
        detray::cesium_iodide_with_ded<scalar>();

    std::string path;
    // Create output directory
    if (output_directory.empty()) {
        path = "";
    } else {
        std::filesystem::create_directories(output_directory);
        path = output_directory + "/";
    }

    // Output Csv file
    std::ofstream helix_rect_file;
    std::ofstream const_rect_file;
    // std::ofstream inhom_rect_no_gradient_file;
    std::ofstream inhom_rect_file;
    std::ofstream inhom_rect_material_file;
    std::ofstream helix_wire_file;
    std::ofstream const_wire_file;
    // std::ofstream inhom_wire_no_gradient_file;
    std::ofstream inhom_wire_file;
    std::ofstream inhom_wire_material_file;
    helix_rect_file.open(path + "helix_rect.csv");
    const_rect_file.open(path + "const_rect.csv");
    // inhom_rect_no_gradient_file.open(path + "inhom_rect_no_gradient.csv");
    inhom_rect_file.open(path + "inhom_rect.csv");
    inhom_rect_material_file.open(path + "inhom_rect_material.csv");
    helix_wire_file.open(path + "helix_wire.csv");
    const_wire_file.open(path + "const_wire.csv");
    // inhom_wire_no_gradient_file.open(path + "inhom_wire_no_gradient.csv");
    inhom_wire_file.open(path + "inhom_wire.csv");
    inhom_wire_material_file.open(path + "inhom_wire_material.csv");

    setup_csv_header_jacobian(helix_rect_file);
    setup_csv_header_jacobian(const_rect_file);
    // setup_csv_header_jacobian(inhom_rect_no_gradient_file);
    setup_csv_header_jacobian(inhom_rect_file);
    setup_csv_header_jacobian(inhom_rect_material_file);
    setup_csv_header_jacobian(helix_wire_file);
    setup_csv_header_jacobian(const_wire_file);
    // setup_csv_header_jacobian(inhom_wire_no_gradient_file);
    setup_csv_header_jacobian(inhom_wire_file);
    setup_csv_header_jacobian(inhom_wire_material_file);

    std::ofstream rect_cov_transport_file;
    rect_cov_transport_file.open(path + "rect_cov_transport.csv");
    setup_csv_header_covariance(rect_cov_transport_file);

    std::ofstream wire_cov_transport_file;
    wire_cov_transport_file.open(path + "wire_cov_transport.csv");
    setup_csv_header_covariance(wire_cov_transport_file);

    // Output Csv file (RK tolerance iteration mode)
    std::vector<std::ofstream> rect_files(log10_tols.size());
    std::vector<std::ofstream> wire_files(log10_tols.size());

    if (rk_tolerance_iterate_mode) {
        for (std::size_t i = 0u; i < log10_tols.size(); i++) {
            const std::string rect_name = "inhom_rect_material_" +
                                          std::to_string(int(log10_tols[i])) +
                                          ".csv";
            const std::string wire_name = "inhom_wire_material_" +
                                          std::to_string(int(log10_tols[i])) +
                                          ".csv";
            rect_files[i].open(path + rect_name);
            wire_files[i].open(path + wire_name);

            setup_csv_header_jacobian(rect_files[i]);
            setup_csv_header_jacobian(wire_files[i]);
        }
    }

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    // Detector types
    using rectangle_telescope = detector<telescope_metadata<rect_type>>;
    using wire_telescope = detector<telescope_metadata<wire_type>>;
    using track_type = free_track_parameters<transform3_type>;

    // Constant magnetic field type
    using const_bfield_t = bfield::const_field_t;

    // Magnetic field map using nearest neightbor interpolation
    using inhom_bfield_t = bfield::inhom_field_t;

    const const_bfield_t const_bfield = bfield::create_const_field(B_z);
    const inhom_bfield_t inhom_bfield = bfield::create_inhom_field();

    // Actor chain type
    using actor_chain_t =
        actor_chain<dtuple, parameter_transporter<transform3_type>,
                    bound_getter<transform3_type>,
                    parameter_resetter<transform3_type>>;

    // Iterate over reference (pilot) tracks for a rectangular telescope
    // geometry and Jacobian calculation
    using uniform_gen_t =
        random_numbers<scalar, std::uniform_real_distribution<scalar>,
                       std::seed_seq>;
    using trk_generator_t = random_track_generator<track_type, uniform_gen_t>;
    trk_generator_t::configuration trk_gen_cfg{};
    trk_gen_cfg.seed(mc_seed);
    trk_gen_cfg.n_tracks(n_tracks + n_skips);
    trk_gen_cfg.phi_range(-constant<scalar>::pi, constant<scalar>::pi);
    trk_gen_cfg.theta_range(0.f, constant<scalar>::pi);
    trk_gen_cfg.mom_range(0.5f * unit<scalar>::GeV, 100.f * unit<scalar>::GeV);
    trk_gen_cfg.origin({0.f, 0.f, 0.f});
    trk_gen_cfg.origin_stddev({0.f * unit<scalar>::mm, 0.f * unit<scalar>::mm,
                               0.f * unit<scalar>::mm});

    // Vectors for dqopdqop relative difference
    std::vector<std::vector<scalar>> dqopdqop_rel_diffs_rect(log10_tols.size());
    std::vector<std::vector<scalar>> dqopdqop_rel_diffs_wire(log10_tols.size());

    bool do_inspect = false;
    if (verbose_lvl >= 4) {
        do_inspect = true;
    }

    // Navigator types
    using rect_navigator_t = navigator<rectangle_telescope>;
    using wire_navigator_t = navigator<wire_telescope>;

    // Stepper types
    using const_field_stepper_t =
        rk_stepper<const_bfield_t::view_t, transform3_type, constrained_step<>,
                   stepper_default_policy>;
    using inhom_field_stepper_t =
        rk_stepper<inhom_bfield_t::view_t, transform3_type, constrained_step<>,
                   stepper_default_policy>;

    // Make four propagators for each case
    using const_field_rect_propagator_t =
        propagator<const_field_stepper_t, rect_navigator_t, actor_chain_t>;
    using inhom_field_rect_propagator_t =
        propagator<inhom_field_stepper_t, rect_navigator_t, actor_chain_t>;

    using const_field_wire_propagator_t =
        propagator<const_field_stepper_t, wire_navigator_t, actor_chain_t>;
    using inhom_field_wire_propagator_t =
        propagator<inhom_field_stepper_t, wire_navigator_t, actor_chain_t>;

    unsigned int track_count = 0u;

    for (const auto track : trk_generator_t{trk_gen_cfg}) {
        mt2.seed(track_count);

        // Pilot track
        detail::helix<transform3_type> helix_bz(track, &B_z);

        // Make a telescope geometry with rectagular surface
        const scalar detector_length = rand_length(mt1);
        const scalar constraint_step_size = detector_length * 1.25f;

        mask<rect_type> rect{0u, detector_length * 1.2f,
                             detector_length * 1.2f};
        mask<wire_type> wire{0u, detector_length * 1.2f,
                             detector_length * 1.2f};

        // Adjust overstep tolerance
        scalar overstep_tol =
            vm["overstep-tolerance-mm"].as<scalar>() * unit<scalar>::mm;
        overstep_tol =
            -math::min(math::abs(overstep_tol), detector_length * 0.5f);

        tel_det_config<rect_type, detail::helix<transform3_type>> rectangle_cfg{
            rect, helix_bz};
        rectangle_cfg.m_envelope = 1000.f * unit<scalar>::mm;
        rectangle_cfg.n_surfaces(2u).length(detector_length);

        auto [rect_det, rect_names] =
            create_telescope_detector(host_mr, rectangle_cfg);
        const auto [euler_rect_initial, shift_rect_initial] =
            tilt_surface<decltype(rect_det),
                         decltype(rect_det)::masks::id::e_rectangle2>(
                rect_det, 0u, helix_bz.dir(0.f));
        const auto [euler_rect_final, shift_rect_final] =
            tilt_surface<decltype(rect_det),
                         decltype(rect_det)::masks::id::e_rectangle2>(
                rect_det, 1u, helix_bz.dir(detector_length));

        // Make a telescope geometry with wire surface
        tel_det_config<wire_type, detail::helix<transform3_type>> wire_cfg{
            wire, helix_bz};
        wire_cfg.m_envelope = 1000.f * unit<scalar>::mm;
        wire_cfg.n_surfaces(2u).length(detector_length);

        auto [wire_det, wire_names] =
            create_telescope_detector(host_mr, wire_cfg);
        const auto [euler_wire_initial, shift_wire_initial] =
            tilt_surface<decltype(wire_det),
                         decltype(wire_det)::masks::id::e_cell_wire>(
                wire_det, 0u, helix_bz.dir(0.f));
        const auto [euler_wire_final, shift_wire_final] =
            tilt_surface<decltype(wire_det),
                         decltype(wire_det)::masks::id::e_cell_wire>(
                wire_det, 1u, helix_bz.dir(detector_length));

        // This IF block should locate after `tilt_surface()` calls for
        // debugging purpose
        if (track_count + 1 <= n_skips) {
            track_count++;
            continue;
        }

        if (verbose_lvl >= 1) {
            track_count++;
            std::cout << "[Event Property]" << std::endl;
            std::cout << "Track ID: " << track_count
                      << "  Number of processed tracks per thread: "
                      << track_count - n_skips << std::endl;
        }
        if (verbose_lvl >= 2) {
            std::cout << "[Detector Property]" << std::endl;
            std::cout << "Path length for the final surface: "
                      << detector_length << std::endl;
            std::cout << "Rect initial surface rotation: ("
                      << euler_rect_initial.alpha << " "
                      << euler_rect_initial.beta << " "
                      << euler_rect_initial.gamma << ")" << std::endl;
            std::cout << "Rect initial surface shift: ("
                      << shift_rect_initial[0u] << " " << shift_rect_initial[1u]
                      << " " << shift_rect_initial[2u] << ")" << std::endl;
            std::cout << "Rect final surface rotation: ("
                      << euler_rect_final.alpha << " " << euler_rect_final.beta
                      << " " << euler_rect_final.gamma << ")" << std::endl;
            std::cout << "Rect final surface shift: (" << shift_rect_final[0u]
                      << " " << shift_rect_final[1u] << " "
                      << shift_rect_final[2u] << ")" << std::endl;
            std::cout << "Wire initial surface rotation: ("
                      << euler_wire_initial.alpha << " "
                      << euler_wire_initial.beta << " "
                      << euler_wire_initial.gamma << ")" << std::endl;
            std::cout << "Wire initial surface shift: ("
                      << shift_wire_initial[0u] << " " << shift_wire_initial[1u]
                      << " " << shift_wire_initial[2u] << ")" << std::endl;
            std::cout << "Wire final surface rotation: ("
                      << euler_wire_final.alpha << " " << euler_wire_final.beta
                      << " " << euler_wire_final.gamma << ")" << std::endl;
            std::cout << "Wire final surface shift: (" << shift_wire_final[0u]
                      << " " << shift_wire_final[1u] << " "
                      << shift_wire_final[2u] << ")" << std::endl;
        }
        if (verbose_lvl >= 3) {
            std::cout << "[Track Property]" << std::endl;
            std::cout << "Phi: " << getter::phi(track.dir()) << std::endl;
            std::cout << "Theta: " << getter::theta(track.dir()) << std::endl;
            std::cout << "Mom: " << track.p() << std::endl;
        }

        /**********************************
         * Rectangluar telescope geometry
         **********************************/

        if (verbose_lvl >= 3 && !skip_rect) {
            std::cout << "Simulating rectangular telescope..." << std::endl;
        }

        // Get initial parameter
        const auto rect_bparam =
            get_initial_parameter<decltype(rect_det),
                                  decltype(rect_det)::masks::id::e_rectangle2>(
                rect_det, track, B_z, helix_tol);

        if (!skip_rect) {

            scalar ref_rel_diff;

            if (rk_tolerance_iterate_mode) {
                for (std::size_t i = 0u; i < log10_tols.size(); i++) {

                    // Rectangle Inhomogeneous field with Material
                    evaluate_jacobian_difference<inhom_field_rect_propagator_t>(
                        track_count, rect_det, detector_length, rect_bparam,
                        inhom_bfield, volume_mat, overstep_tol, on_surface_tol,
                        std::pow(10.f, log10_tols[i]), constraint_step_size,
                        h_sizes_rect, rect_files[i], ref_rel_diff, true,
                        do_inspect);

                    dqopdqop_rel_diffs_rect[i].push_back(ref_rel_diff);
                }
            } else if (!rk_tolerance_iterate_mode) {

                // For helix
                evaluate_jacobian_difference_helix<
                    decltype(rect_det),
                    decltype(rect_det)::masks::id::e_rectangle2>(
                    track_count, rect_det, rect_bparam, B_z, h_sizes_rect,
                    helix_rect_file, helix_tol);

                // Rect Const field
                evaluate_jacobian_difference<const_field_rect_propagator_t>(
                    track_count, rect_det, detector_length, rect_bparam,
                    const_bfield, vacuum<scalar>(), overstep_tol,
                    on_surface_tol, rk_tol, constraint_step_size, h_sizes_rect,
                    const_rect_file, ref_rel_diff);

                /*
                // Rect Inhomogeneous field with no gradient
                evaluate_jacobian_difference<inhom_field_rect_propagator_t>(
                    track_count, rect_det, rect_bparam, inhom_bfield,
                    vacuum<scalar>(), overstep_tol, on_surface_tol, rk_tol,
                    constraint_step_size, h_sizes_rect,
                    inhom_rect_no_gradient_file, ref_rel_diff, false);
                */

                // Rect Inhomogeneous field
                evaluate_jacobian_difference<inhom_field_rect_propagator_t>(
                    track_count, rect_det, detector_length, rect_bparam,
                    inhom_bfield, vacuum<scalar>(), overstep_tol,
                    on_surface_tol, rk_tol, constraint_step_size, h_sizes_rect,
                    inhom_rect_file, ref_rel_diff);

                // Rectangle Inhomogeneous field with Material
                evaluate_jacobian_difference<inhom_field_rect_propagator_t>(
                    track_count, rect_det, detector_length, rect_bparam,
                    inhom_bfield, volume_mat, overstep_tol, on_surface_tol,
                    rk_tol, constraint_step_size, h_sizes_rect,
                    inhom_rect_material_file, ref_rel_diff);

                // Rectangle Inhomogeneous field with Material (Covariance
                // transport)
                evaluate_covariance_transport<inhom_field_rect_propagator_t>(
                    track_count, rect_det, detector_length, rect_bparam,
                    inhom_bfield, volume_mat, overstep_tol, on_surface_tol,
                    rk_tol, constraint_step_size, rect_cov_transport_file);
            }
        }

        /**********************************
         * Wire telescope geometry
         **********************************/

        if (verbose_lvl >= 3 && !skip_wire) {
            std::cout << "Simulating wire telescope..." << std::endl;
        }

        // Get initial parameter
        const auto wire_bparam =
            get_initial_parameter<decltype(wire_det),
                                  decltype(wire_det)::masks::id::e_cell_wire>(
                wire_det, track, B_z, helix_tol);

        if (!skip_wire) {

            scalar ref_rel_diff;

            if (rk_tolerance_iterate_mode) {
                for (std::size_t i = 0u; i < log10_tols.size(); i++) {
                    // Wire Inhomogeneous field with Material
                    evaluate_jacobian_difference<inhom_field_wire_propagator_t>(
                        track_count, wire_det, detector_length, wire_bparam,
                        inhom_bfield, volume_mat, overstep_tol, on_surface_tol,
                        std::pow(10.f, log10_tols[i]), constraint_step_size,
                        h_sizes_wire, wire_files[i], ref_rel_diff, true,
                        do_inspect);

                    dqopdqop_rel_diffs_wire[i].push_back(ref_rel_diff);
                }
            } else if (!rk_tolerance_iterate_mode) {

                // For helix
                evaluate_jacobian_difference_helix<
                    decltype(wire_det),
                    decltype(wire_det)::masks::id::e_cell_wire>(
                    track_count, wire_det, wire_bparam, B_z, h_sizes_wire,
                    helix_wire_file, helix_tol);

                // Wire Const field
                evaluate_jacobian_difference<const_field_wire_propagator_t>(
                    track_count, wire_det, detector_length, wire_bparam,
                    const_bfield, vacuum<scalar>(), overstep_tol,
                    on_surface_tol, rk_tol, constraint_step_size, h_sizes_wire,
                    const_wire_file, ref_rel_diff);

                /*
                // Wire Inhomogeneous field with no gradient
                evaluate_jacobian_difference<inhom_field_wire_propagator_t>(
                    track_count, wire_det, wire_bparam, inhom_bfield,
                    vacuum<scalar>(), overstep_tol, on_surface_tol, rk_tol,
                    constraint_step_size, h_sizes_wire,
                    inhom_wire_no_gradient_file, ref_rel_diff, false);
                */

                // Wire Inhomogeneous field
                evaluate_jacobian_difference<inhom_field_wire_propagator_t>(
                    track_count, wire_det, detector_length, wire_bparam,
                    inhom_bfield, vacuum<scalar>(), overstep_tol,
                    on_surface_tol, rk_tol, constraint_step_size, h_sizes_wire,
                    inhom_wire_file, ref_rel_diff);

                // Wire Inhomogeneous field with Material
                evaluate_jacobian_difference<inhom_field_wire_propagator_t>(
                    track_count, wire_det, detector_length, wire_bparam,
                    inhom_bfield, volume_mat, overstep_tol, on_surface_tol,
                    rk_tol, constraint_step_size, h_sizes_wire,
                    inhom_wire_material_file, ref_rel_diff);

                // Wire Inhomogeneous field with Material (Covariance transport)
                evaluate_covariance_transport<inhom_field_wire_propagator_t>(
                    track_count, wire_det, detector_length, wire_bparam,
                    inhom_bfield, volume_mat, overstep_tol, on_surface_tol,
                    rk_tol, constraint_step_size, wire_cov_transport_file);
            }
        }
    }

    if (rk_tolerance_iterate_mode) {
        for (std::size_t i = 0u; i < log10_tols.size(); i++) {
            EXPECT_EQ(dqopdqop_rel_diffs_rect[i].size(), n_tracks);
            EXPECT_EQ(dqopdqop_rel_diffs_wire[i].size(), n_tracks);

            EXPECT_GE(statistics::mean(dqopdqop_rel_diffs_rect[i]), 1e-10f);
            EXPECT_LE(statistics::mean(dqopdqop_rel_diffs_rect[i]), 1e-2f);
            EXPECT_GE(statistics::mean(dqopdqop_rel_diffs_wire[i]), 1e-10f);
            EXPECT_LE(statistics::mean(dqopdqop_rel_diffs_wire[i]), 1e-2f);
        }
    }

    // Close files
    helix_rect_file.close();
    const_rect_file.close();
    // inhom_rect_no_gradient_file.close();
    inhom_rect_file.close();
    inhom_rect_material_file.close();

    helix_wire_file.close();
    const_wire_file.close();
    // inhom_wire_no_gradient_file.close();
    inhom_wire_file.close();
    inhom_wire_material_file.close();

    rect_cov_transport_file.close();
    wire_cov_transport_file.close();

    if (rk_tolerance_iterate_mode) {
        for (std::size_t i = 0u; i < log10_tols.size(); i++) {
            rect_files[i].close();
            wire_files[i].close();
        }
    }
}
