/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "tests/common/benchmark/benchmark_propagation.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace detray;

using scalar_t = detray::scalar;

int main(int argc, char** argv) {

    detray::benchmark_base::configuration cfg{};
    cfg.do_sleep(false);

    // Prepare tracks
    const unsigned int phi_steps{5u};
    const unsigned int theta_steps{5u};
    const unsigned int n_trks{phi_steps * theta_steps};

    // Add additional tracks for warmup
    const unsigned int add_trks{
        cfg.do_warmup()
            ? static_cast<unsigned int>(std::ceil(
                  -0.5f * static_cast<scalar_t>(phi_steps + theta_steps) +
                  std::sqrt(static_cast<scalar_t>(cfg.n_warmup()) +
                            0.25f * std::pow(phi_steps + theta_steps, 2))))
            : 0u};

    // Detector
    toy_detector_fixture toy_det_fixture{};

    // Generate tracks
    using track_generator_t = uniform_track_generator<
        free_track_parameters<toy_detector_fixture::transform3_type>>;
    track_generator_t::configuration trk_cfg{};
    trk_cfg.theta_steps(phi_steps + add_trks)
        .phi_steps(theta_steps + add_trks)
        .p_mag(10.f * unit<scalar_t>::GeV);
    std::vector<free_track_parameters<toy_detector_fixture::transform3_type>>
        tracks;

    cfg.n_samples(n_trks).n_warmup(
        std::ceil(0.1f * static_cast<scalar_t>(n_trks)));

    fill_tracks<track_generator_t>(tracks, trk_cfg);

    //
    // Register all benchmarks
    //

    // Default navigation policy
    rkn_propagation_bm<always_init> default_rkn_toygeo{cfg};
    default_rkn_toygeo.config()
        .detector_fixture()
        .config()
        .n_barrel_layers(0u)
        .n_endcap_layers(0u);

    // Trivial geometry
    std::string name{default_rkn_toygeo.name() + "_only_beampipe"};

    ::benchmark::RegisterBenchmark(name.c_str(), default_rkn_toygeo, tracks)
        ->MeasureProcessCPUTime();

    // Varying number of barrel layers, no endcaps
    for (unsigned int n{1u}; n < 5u; ++n) {
        name = default_rkn_toygeo.name() + "_" + std::to_string(n) + "brl_0edc";
        default_rkn_toygeo.config()
            .detector_fixture()
            .config()
            .n_barrel_layers(n)
            .n_endcap_layers(0u);

        ::benchmark::RegisterBenchmark(name.c_str(), default_rkn_toygeo, tracks)
            ->MeasureProcessCPUTime();
    }
    // Varying number of endcap layers
    for (unsigned int n{1u}; n < 8u; ++n) {
        name = default_rkn_toygeo.name() + "_4brl_" + std::to_string(n) + "edc";
        default_rkn_toygeo.config()
            .detector_fixture()
            .config()
            .n_barrel_layers(4u)
            .n_endcap_layers(n);

        ::benchmark::RegisterBenchmark(name.c_str(), default_rkn_toygeo, tracks)
            ->MeasureProcessCPUTime();
    }

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
}