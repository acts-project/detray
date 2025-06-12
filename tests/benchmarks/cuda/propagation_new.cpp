/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/device/cuda/propagator.hpp"
#include "detray/benchmarks/propagation_benchmark_utils.hpp"
#include "detray/benchmarks/types.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/track_generators.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <chrono>
#include <ctime>
#include <iostream>
#include <ratio>
#include <string>

using namespace detray;

int main(int argc, char** argv) {

    using metadata_t = benchmarks::toy_metadata;
    using toy_detector_t = detector<metadata_t>;
    using algebra_t = typename toy_detector_t::algebra_type;
    using scalar = dscalar<algebra_t>;
    using vector3 = dvector3D<algebra_t>;

    using free_track_parameters_t = free_track_parameters<algebra_t>;
    using uniform_gen_t =
        detail::random_numbers<scalar, std::uniform_real_distribution<scalar>>;
    using track_generator_t =
        random_track_generator<free_track_parameters_t, uniform_gen_t>;
    using field_bknd_t = bfield::const_bknd_t<benchmarks::scalar>;

    vecmem::host_memory_resource host_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    //
    // Configuration
    //

    // Constant magnetic field
    vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};

    // Configure toy detector
    toy_det_config<scalar> toy_cfg{};
    toy_cfg.use_material_maps(false).n_brl_layers(4u).n_edc_layers(7u);

    std::cout << toy_cfg << std::endl;

    // Configure propagation
    propagation::config prop_cfg{};
    prop_cfg.navigation.search_window = {3u, 3u};

    std::cout << prop_cfg << std::endl;

    //
    // Prepare data
    //
    // Generate track sample for strong scaling
    track_generator_t::configuration trk_cfg{};
    trk_cfg.n_tracks(10u);
    trk_cfg.seed(detail::random_numbers<scalar>::default_seed());

    track_generator_t trk_gen{trk_cfg};

    dvector<free_track_parameters_t> single_sample =
        detray::benchmarks::generate_tracks(&host_mr, trk_gen);

    const auto [toy_det, names] =
        build_toy_detector<algebra_t>(host_mr, toy_cfg);

    auto bfield = create_const_field<scalar>(B);

    pointwise_material_interactor<algebra_t>::state interactor_state{};

    auto actor_states = detail::make_tuple<dtuple>(interactor_state);

    //
    // Register benchmarks
    //
    std::cout << "----------------------\n"
              << "Propagation Test\n"
              << "----------------------\n\n";

    using navigator_t = navigator_type<metadata_t>;
    using stepper_t = stepper_type<metadata_t, field_bknd_t>;
    using actor_chain_t = default_chain<algebra_t>;

    prop_cfg.stepping.do_covariance_transport = true;
    cuda_propagation<navigator_t, stepper_t, actor_chain_t> propagator{
        prop_cfg};

    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();
    propagator(&dev_mr, &toy_det, &bfield, &single_sample, &actor_states);
    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();

    const auto time_span =
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    std::cout << "It took me " << time_span << " for " << trk_cfg.n_tracks()
              << " tracks" << std::endl;
}
