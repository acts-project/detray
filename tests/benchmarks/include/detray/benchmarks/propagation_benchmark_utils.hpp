/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/containers.hpp"
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
#include <utility>
#include <vector>

namespace detray {

/// Which propagate function to run
enum class propagate_option {
    e_unsync = 0,
    e_sync = 1,
};

/// Define propagator type
template <typename stepper_t, typename detector_t,
          typename actor_chain_t = actor_chain<>>
using propagator_t =
    propagator<stepper_t, navigator<detector_t>, actor_chain_t>;

/// Precompute the tracks
template <typename track_generator_t>
inline auto generate_tracks(
    vecmem::memory_resource *mr,
    const typename track_generator_t::configuration &cfg = {}) {

    // Track collection
    dvector<typename track_generator_t::track_type> tracks(mr);

    // Iterate through uniformly distributed momentum directions
    for (auto track : track_generator_t{cfg}) {
        // Put it into vector of trajectories
        tracks.push_back(track);
    }

    return tracks;
}

/// Tie the actor states for the propagation
template <typename propagator_t, std::size_t... I>
inline constexpr auto setup_actor_states(
    typename propagator_t::actor_chain_type::state_tuple &input_actor_states,
    std::index_sequence<I...>) {

    // Empty actor chain
    if constexpr (sizeof...(I) == 0u) {
        return typename propagator_t::actor_chain_type::state{};
    } else {
        return detray::tie(detail::get<I>(input_actor_states)...);
    }
}

/// Register a propagation benchmark of type @tparam benchmark_t
///
/// @tparam benchmark_t the propagation benchmark functor
/// @tparam stepper_t the stepper to use fro track parameter transport
/// @tparam actor_chain_t types of actors
/// @tparam detector_t detector type
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
template <template <typename, typename> class benchmark_t, typename stepper_t,
          typename actor_chain_t, typename detector_t, typename bfield_t>
inline void register_benchmark(
    const std::string &name, benchmark_base::configuration &bench_cfg,
    propagation::config &prop_cfg, const detector_t &det,
    const bfield_t &bfield, typename actor_chain_t::state_tuple &actor_states,
    dvector<free_track_parameters<typename detector_t::algebra_type>> &tracks,
    const std::vector<int> &n_samples = {10000}) {

    for (int n : n_samples) {

        bench_cfg.n_samples(n);

        // Without covariance transport
        benchmark_t<propagator_t<stepper_t, detector_t, actor_chain_t>,
                    bfield_t>
            prop_benchmark{bench_cfg};
        prop_benchmark.config().propagation() = prop_cfg;

        std::string bench_name = prop_benchmark.config().name() + "_" + name +
                                 "_" + std::to_string(n) + "_TRACKS";

        std::cout << bench_name << "\n" << bench_cfg;

        ::benchmark::RegisterBenchmark(bench_name.c_str(), prop_benchmark,
                                       &tracks, &det, &bfield, &actor_states);
        //->MeasureProcessCPUTime();
    }
}

}  // namespace detray
