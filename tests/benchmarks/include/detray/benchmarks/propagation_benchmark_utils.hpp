/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/tuple.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <type_traits>
#include <utility>
#include <vector>

namespace detray::benchmarks {

/// Which propagate function to run
enum class propagation_opt {
    e_unsync = 0,
    e_sync = 1,
};

/// @returns the default track generation configuration for detray benchmarks
template <typename track_generator_t>
inline typename track_generator_t::configuration get_default_trk_gen_config(
    const std::vector<int> &n_tracks) {

    using track_t = typename track_generator_t::track_type;
    using scalar_t = dscalar<typename track_t::algebra_type>;

    int n_trks{*std::ranges::max_element(n_tracks)};

    // Generate tracks
    typename track_generator_t::configuration trk_cfg{};
    trk_cfg.n_tracks(static_cast<std::size_t>(n_trks));
    trk_cfg.randomize_charge(true);
    trk_cfg.phi_range(-constant<scalar_t>::pi, constant<scalar_t>::pi);
    trk_cfg.eta_range(-3.f, 3.f);
    trk_cfg.mom_range(1.f * unit<scalar_t>::GeV, 100.f * unit<scalar_t>::GeV);
    trk_cfg.origin(0.f, 0.f, 0.f);
    trk_cfg.origin_stddev(0.f, 0.f, 0.f);

    return trk_cfg;
}

/// Precompute the tracks
///
/// @param mr memory resource to allocate the track vector
/// @param cfg the configuration of the track generator
/// @param do_sort sort the tracks by theta angle
template <typename track_generator_t>
inline auto generate_tracks(
    vecmem::memory_resource *mr,
    const typename track_generator_t::configuration &cfg = {},
    bool do_sort = true) {

    using track_t = typename track_generator_t::track_type;
    using scalar_t = dscalar<typename track_t::algebra_type>;

    // Track collection
    dvector<track_t> tracks(mr);

    // Iterate through uniformly distributed momentum directions
    for (auto track : track_generator_t{cfg}) {
        // Put it into vector of trajectories
        tracks.push_back(track);
    }

    if (do_sort) {
        // Sort by theta angle
        const auto traj_comp = [](const auto &lhs, const auto &rhs) {
            constexpr auto pi_2{constant<scalar_t>::pi_2};
            return math::fabs(pi_2 - vector::theta(lhs.dir())) <
                   math::fabs(pi_2 - vector::theta(rhs.dir()));
        };

        std::ranges::sort(tracks, traj_comp);
    }

    return tracks;
}

/// Generate as many samples of track states as there are entries in the
/// @param n_tracks vector.
template <typename track_generator_t>
inline auto generate_track_samples(
    vecmem::memory_resource *mr, const std::vector<int> &n_tracks,
    typename track_generator_t::configuration &cfg = {}, bool do_sort = true) {

    using track_t = typename track_generator_t::track_type;

    std::vector<dvector<track_t>> track_samples{};
    track_samples.reserve(n_tracks.size());

    auto tmp_cfg{cfg};
    for (const int n : n_tracks) {
        tmp_cfg.n_tracks(static_cast<std::size_t>(n));
        track_samples.push_back(
            generate_tracks<track_generator_t>(mr, tmp_cfg, do_sort));
    }

    return track_samples;
}

/// Register a propagation benchmark of type @tparam benchmark_t
///
/// @tparam benchmark_t the propagation benchmark functor
/// @tparam propagator_t full propagator type
/// @tparam detector_t host detector type
/// @tparam bfield_t covfie magnetic field type
///
/// @param name name for the benchmark
/// @param bench_cfg basic benchmark configuration
/// @param prop_cfg propagation configuration
/// @param det the detector
/// @param bfield the covfie field
/// @param actor_states tuple that contains all actor states (same order as in
///                     actor_chain_t)
/// @param tracks the pre-computed test tracks
/// @param n_samples the number of track to run
template <template <typename, typename, detray::benchmarks::propagation_opt>
          class benchmark_t,
          typename propagator_t, typename detector_t, typename bfield_bknd_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
inline void register_benchmark(
    const std::string &name, benchmark_base::configuration &bench_cfg,
    propagation::config &prop_cfg, const detector_t &det, bfield_bknd_t &bfield,
    typename propagator_t::actor_chain_type::state_tuple *actor_states,
    std::vector<
        dvector<free_track_parameters<typename detector_t::algebra_type>>>
        &track_samples,
    const std::vector<int> &n_samples = {10000},
    vecmem::memory_resource *dev_mr = nullptr) {

    using algebra_t = typename detector_t::algebra_type;
    using propagation_benchmark_t =
        benchmark_t<propagator_t, bfield_bknd_t, kOPT>;

    assert(track_samples.size() == n_samples.size());

    for (const auto [i, n] : detray::views::enumerate(n_samples)) {

        auto &tracks = track_samples[i];
        assert(static_cast<std::size_t>(n) <= tracks.size());

        bench_cfg.n_samples(n);

        typename propagation_benchmark_t::configuration prop_bm_cfg{bench_cfg};
        prop_bm_cfg.propagation() = prop_cfg;

        // Configure the benchmark
        propagation_benchmark_t prop_benchmark{prop_bm_cfg};

        std::string bench_name = prop_benchmark.config().name() + "_" + name +
                                 "_" + std::to_string(n) + "_TRACKS";

        std::cout << bench_name << "\n" << bench_cfg;

        if constexpr (std::is_invocable_v<
                          decltype(prop_benchmark), ::benchmark::State &,
                          dvector<free_track_parameters<algebra_t>> *,
                          const detector_t *, const bfield_bknd_t *,
                          typename propagator_t::actor_chain_type::state_tuple
                              *>) {
            // Cpu benchmark
            ::benchmark::RegisterBenchmark(bench_name.c_str(), prop_benchmark,
                                           &tracks, &det, &bfield,
                                           actor_states);
            //->MeasureProcessCPUTime();
        } else {
            // Device benchmark
            ::benchmark::RegisterBenchmark(bench_name.c_str(), prop_benchmark,
                                           dev_mr, &tracks, &det, &bfield,
                                           actor_states);
            //->MeasureProcessCPUTime();
        }
    }
}

/// Register a propagation benchmark of type @tparam benchmark_t
///
/// @tparam benchmark_t the propagation benchmark functor
/// @tparam stepper_t the stepper to use fro track parameter transport
/// @tparam actor_chain_t types of actors
template <template <typename, typename, detray::benchmarks::propagation_opt>
          class benchmark_t,
          typename stepper_t, typename actor_chain_t, typename detector_t,
          typename bfield_bknd_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
inline void register_benchmark(
    const std::string &name, benchmark_base::configuration &bench_cfg,
    propagation::config &prop_cfg, const detector_t &det, bfield_bknd_t &bfield,
    typename actor_chain_t::state_tuple *actor_states,
    std::vector<
        dvector<free_track_parameters<typename detector_t::algebra_type>>>
        &tracks,
    const std::vector<int> &n_samples = {10000}) {

    using propagator_t =
        propagator<stepper_t, navigator<detector_t>, actor_chain_t>;
    register_benchmark<benchmark_t, propagator_t, detector_t, bfield_bknd_t,
                       kOPT>(name, bench_cfg, prop_cfg, det, bfield,
                             actor_states, tracks, n_samples, nullptr);
}

}  // namespace detray::benchmarks
