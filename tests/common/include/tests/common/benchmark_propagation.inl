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

using track_generator_t = uniform_track_generator<free_track_parameters<transform3>>;

int main() {

    detray::benchmark::config cfg{};
    cfg.do_warmup(false).do_sleep(false);

    // Prepare tracks
    const unsigned int phi_steps{5u};
    const unsigned int theta_steps{5u};
    const unsigned int n_trks{phi_steps * theta_steps};
    const unsigned int add_trks{cfg.do_warmup() ? static_cast<unsigned int>(std::ceil(-0.5f * (phi_steps + theta_steps) + std::sqrt(cfg.n_warmup() + 0.25f * std::pow(phi_steps + theta_steps, 2)))) : 0u};

    track_generator_t::configuration trk_cfg{};
    trk_cfg.theta_steps(phi_steps + add_trks).phi_steps(theta_steps + add_trks).p_mag(10.f * unit<scalar>::GeV);
    std::vector<free_track_parameters<transform3>> tracks;

    fill_tracks<track_generator_t>(tracks, trk_cfg);

    //
    // Register all benchmarks
    //

    // Default navigation policy
    rkn_toy_bm<always_init> default_rkn_toygeo{cfg};

    // Trivial geometry
    std::string name{rkn_toy_bm<>::name + "_only_beampipe"};
    default_rkn_toygeo.config().n_barrel_layers(0u).n_endcap_layers(0u);

    ::benchmark::RegisterBenchmark(name.c_str(), default_rkn_toygeo, tracks, n_trks)
        ->Repetitions(5)
        ->DisplayAggregatesOnly(true);

    // Varying number of barrel layers, no endcaps
    for (unsigned int n{1u}; n < 5u; ++n) {
        name = rkn_toy_bm<>::name + "_" + std::to_string(n) + "brl_0edc";
        default_rkn_toygeo.config().n_barrel_layers(n).n_endcap_layers(0u);

        ::benchmark::RegisterBenchmark(name.c_str(), default_rkn_toygeo, tracks, n_trks)
        ->Repetitions(5)
        ->DisplayAggregatesOnly(true);
    }
    // Varying number of endcap layers
    for (unsigned int n{1u}; n < 8u; ++n) {
        name = rkn_toy_bm<>::name + "_4brl_" + std::to_string(n) + "edc";
        default_rkn_toygeo.config().n_barrel_layers(4u).n_endcap_layers(n);

        ::benchmark::RegisterBenchmark(name.c_str(), default_rkn_toygeo, tracks,  n_trks)
        ->Repetitions(5)
        ->DisplayAggregatesOnly(true);
    }

    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
}