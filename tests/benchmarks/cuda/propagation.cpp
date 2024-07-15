/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/benchmarks/device/cuda/propagation_benchmark.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/detectors/create_wire_chamber.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <string>

using namespace detray;

int main(int argc, char** argv) {

    using toy_detector_t = detector<toy_metadata>;
    using algebra_t = typename toy_detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;
    using free_track_parameters_t = free_track_parameters<algebra_t>;
    using uniform_gen_t =
        detail::random_numbers<scalar_t,
                               std::uniform_real_distribution<scalar_t>>;
    using track_generator_t =
        random_track_generator<free_track_parameters_t, uniform_gen_t>;
    using field_bknd_t = bfield::const_bknd_t;

    vecmem::host_memory_resource host_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    //
    // Configuration
    //

    // Constant magnetic field
    vector3_t B{0.f, 0.f, 2.f * unit<scalar_t>::T};

    // Configure toy detector
    toy_det_config toy_cfg{};
    toy_cfg.use_material_maps(true).n_brl_layers(4u).n_edc_layers(7u);

    std::cout << toy_cfg << std::endl;

    // Configure wire chamber
    wire_chamber_config wire_chamber_cfg{};
    wire_chamber_cfg.half_z(500.f * unit<scalar>::mm);

    std::cout << wire_chamber_cfg << std::endl;

    // Configure propagation
    propagation::config prop_cfg{};
    prop_cfg.navigation.search_window = {3u, 3u};

    std::cout << prop_cfg << std::endl;

    // Benchmark config
    detray::benchmark_base::configuration bench_cfg{};

    std::vector<int> n_tracks{8 * 8,   16 * 16,   32 * 32,
                              64 * 64, 128 * 128, 256 * 256};

    int n_trks{*std::max_element(std::begin(n_tracks), std::end(n_tracks))};
    std::cout << n_trks << std::endl;
    bench_cfg.n_warmup(
        static_cast<int>(std::ceil(0.1f * static_cast<float>(n_trks))));
    // Add tracks for warmup
    n_trks += bench_cfg.do_warmup() ? bench_cfg.n_warmup() : 0;

    // Generate tracks
    track_generator_t::configuration trk_cfg{};
    trk_cfg.n_tracks(n_trks)
        .randomize_charge(true)
        .eta_range(-4.f, 4.f)
        .pT_range(0.5f * unit<scalar_t>::GeV, 100.f * unit<scalar_t>::GeV);

    std::cout << trk_cfg << std::endl;

    //
    // Prepare data
    //
    auto tracks = generate_tracks<track_generator_t>(&host_mr, trk_cfg);

    const auto [toy_det, names] = build_toy_detector(host_mr, toy_cfg);
    const auto [wire_chamber, _] =
        create_wire_chamber(host_mr, wire_chamber_cfg);

    auto bfield = bfield::create_const_field(B);

    //
    // Register benchmarks
    //
    std::cout << "Propagation Benchmarks\n"
              << "----------------------\n\n";

    prop_cfg.stepping.do_covariance_transport = true;
    register_benchmark<cuda_propagation_bm,
                       cuda_propagator_type<toy_metadata, field_bknd_t>>(
        "TOY_DETECTOR_W_COV_TRANSPORT", bench_cfg, prop_cfg, toy_det, bfield,
        tracks, n_tracks, &host_mr, &dev_mr);

    prop_cfg.stepping.do_covariance_transport = true;
    register_benchmark<cuda_propagation_bm,
                       cuda_propagator_type<default_metadata, field_bknd_t>>(
        "WIRE_CHAMBER_W_COV_TRANSPORT", bench_cfg, prop_cfg, wire_chamber,
        bfield, tracks, n_tracks, &host_mr, &dev_mr);

    // Run benchmarks
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
}
