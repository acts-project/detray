/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/benchmarks/propagation_benchmark_config.hpp"
#include "detray/benchmarks/propagation_benchmark_utils.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/containers.hpp"
#include "detray/tracks/tracks.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <algorithm>
#include <cassert>
#include <random>
#include <string>

namespace detray {

template <typename propagator_t, typename bfield_t,
          propagate_option opt = propagate_option::e_unsync>
struct propagation_bm : public benchmark_base {
    /// Detector dependent types
    using algebra_t = typename propagator_t::detector_type::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    /// Local configuration type
    using configuration = propagation_benchmark_config;

    /// The benchmark configuration
    configuration m_cfg{};

    /// Default construction
    propagation_bm() = default;

    /// Construct from an externally provided configuration @param cfg
    propagation_bm(configuration cfg) : m_cfg{cfg} {}

    /// @return the benchmark configuration
    configuration &config() { return m_cfg; }

    /// Prepare data and run benchmark loop
    inline void operator()(::benchmark::State &state,
                           dvector<free_track_parameters<algebra_t>> *tracks,
                           const typename propagator_t::detector_type *det,
                           const bfield_t *bfield,
                           typename propagator_t::actor_chain_type::state_tuple
                               *input_actor_states) const {

        using actor_states_t =
            typename propagator_t::actor_chain_type::state_tuple;

        const int n_samples{m_cfg.benchmark().n_samples()};
        const int n_warmup{m_cfg.benchmark().n_warmup()};

        assert(static_cast<std::size_t>(n_samples + n_warmup) <=
               tracks->size());

        // Shuffle the sample
        std::random_device rd;
        std::mt19937 gen(rd());

        std::shuffle(tracks->begin(), tracks->end(), gen);

        // Create propagator
        propagator_t p{m_cfg.propagation()};

        // Warm-up
        if (m_cfg.benchmark().do_warmup()) {
#pragma omp parallel for
            for (int i = 0; i < n_warmup; ++i) {
                // Fresh copy of actor states
                actor_states_t actor_state_tuple(*input_actor_states);
                // Tuple of references to pass to the propagator
                typename propagator_t::actor_chain_type::state actor_states =
                    setup_actor_states<propagator_t>(
                        actor_state_tuple,
                        std::make_integer_sequence<
                            std::size_t,
                            detail::tuple_size_v<actor_states_t>>{});

                typename propagator_t::state p_state((*tracks)[i], *bfield,
                                                     *det);

                // Run propagation
                if constexpr (opt == propagate_option::e_unsync) {
                    ::benchmark::DoNotOptimize(
                        p.propagate(p_state, actor_states));
                } else if constexpr (opt == propagate_option::e_sync) {
                    ::benchmark::DoNotOptimize(
                        p.propagate_sync(p_state, actor_states));
                }
            }
        }

        // Run the benchmark
        for (auto _ : state) {
#pragma omp parallel for
            for (int i = n_warmup; i < n_samples + n_warmup; ++i) {
                // Fresh copy of actor states
                actor_states_t actor_state_tuple(*input_actor_states);
                // Tuple of references to pass to the propagator
                typename propagator_t::actor_chain_type::state actor_states =
                    setup_actor_states<propagator_t>(
                        actor_state_tuple,
                        std::make_integer_sequence<
                            std::size_t,
                            detail::tuple_size_v<actor_states_t>>{});

                typename propagator_t::state p_state((*tracks)[i], *bfield,
                                                     *det);

                // Run propagation
                if constexpr (opt == propagate_option::e_unsync) {
                    ::benchmark::DoNotOptimize(
                        p.propagate(p_state, actor_states));
                } else if constexpr (opt == propagate_option::e_sync) {
                    ::benchmark::DoNotOptimize(
                        p.propagate_sync(p_state, actor_states));
                }
            }
        }
        // Report throughput
        state.counters["TracksPropagated"] = benchmark::Counter(
            static_cast<double>(n_samples), benchmark::Counter::kIsRate);
    }
};

}  // namespace detray
