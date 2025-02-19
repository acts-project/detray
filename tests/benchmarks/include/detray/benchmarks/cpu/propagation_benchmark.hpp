/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/tracks/tracks.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/benchmarks/propagation_benchmark_config.hpp"
#include "detray/benchmarks/propagation_benchmark_utils.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

#ifdef _OPENMP
// openMP include
#include <omp.h>
#endif

// System include(s)
#include <algorithm>
#include <cassert>
#include <ranges>
#include <string>

namespace detray::benchmarks {

template <typename propagator_t, typename bfield_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
struct host_propagation_bm : public benchmark_base {
    /// Detector dependent types
    using algebra_t = typename propagator_t::detector_type::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    /// Local configuration type
    using configuration = propagation_benchmark_config;

    /// The benchmark configuration
    configuration m_cfg{};

    /// Default construction
    host_propagation_bm() = default;

    /// Construct from an externally provided configuration @param cfg
    explicit host_propagation_bm(const configuration &cfg) : m_cfg{cfg} {}

    /// @return the benchmark configuration
    configuration &config() { return m_cfg; }

    /// Prepare data and run benchmark loop
    inline void operator()(
        ::benchmark::State &state,
        const dvector<free_track_parameters<algebra_t>> *tracks,
        const typename propagator_t::detector_type *det, const bfield_t *bfield,
        const typename propagator_t::actor_chain_type::state_tuple
            *input_actor_states,
        [[maybe_unused]] const int n_threads,
        [[maybe_unused]] const int max_chunk_size,
        [[maybe_unused]] const int thread_schedule) const {
        using actor_chain_t = typename propagator_t::actor_chain_type;
        using actor_states_t = typename actor_chain_t::state_tuple;

        assert(tracks != nullptr);
        assert(det != nullptr);
        assert(bfield != nullptr);
        assert(input_actor_states != nullptr);

        const int n_samples{m_cfg.benchmark().n_samples()};
        const int n_warmup{m_cfg.benchmark().n_warmup()};

        assert(static_cast<std::size_t>(n_samples) <= tracks->size());

#ifdef _OPENMP
        // Set the number of threads fo the openMP parallel regions
        omp_set_num_threads(n_threads);
        // Clamp chunk size to [1, max_chunk_size]
        int chunk_size{
            math::min(static_cast<int>(n_samples / n_threads), max_chunk_size)};
        chunk_size = math::max(chunk_size, 1);
        omp_set_schedule(static_cast<omp_sched_t>(thread_schedule), chunk_size);
#ifndef NDEBUG
        std::cout << "No. tracks " << n_samples << std::endl;
        std::cout << "No. threads " << n_threads << std::endl;
        std::cout << "Schedule type " << thread_schedule << std::endl;
        std::cout << "Chunk size " << chunk_size << std::endl;
#endif
#endif

        // Create propagator
        propagator_t p{m_cfg.propagation()};

        // Call the host propagation
        auto run_propagation = [&p, det, bfield, input_actor_states](
                                   const free_track_parameters<algebra_t>
                                       &track) {
            // Fresh copy of actor states
            actor_states_t actor_states(*input_actor_states);
            // Tuple of references to pass to the propagator
            typename actor_chain_t::state_ref_tuple actor_state_refs =
                actor_chain_t::setup_actor_states(actor_states);

            typename propagator_t::state p_state(track, *bfield, *det);
            // Particle hypothesis
            auto &ptc = p_state._stepping.particle_hypothesis();
            p_state.set_particle(update_particle_hypothesis(ptc, track));

            // Run propagation
            if constexpr (kOPT ==
                          detray::benchmarks::propagation_opt::e_unsync) {
                ::benchmark::DoNotOptimize(
                    p.propagate(p_state, actor_state_refs));
            } else if constexpr (kOPT ==
                                 detray::benchmarks::propagation_opt::e_sync) {
                ::benchmark::DoNotOptimize(
                    p.propagate_sync(p_state, actor_state_refs));
            }
            assert(p.propagate_is_complete(p_state));
        };

        // Warm-up
        if (m_cfg.benchmark().do_warmup()) {
            assert(n_warmup > 0);
            int stride{n_samples / n_warmup};
            stride = (stride == 0) ? 10 : stride;
            assert(stride > 0);

#pragma omp parallel for
            for (int i = 0; i < n_samples; i += stride) {
                // The track gets copied into the stepper state, so that the
                // original track sample vector remains unchanged
                run_propagation((*tracks)[static_cast<std::size_t>(i)]);
            }
        } else {
            std::cout << "WARNING: Running host benchmarks without warmup"
                      << std::endl;
        }

        // Run the benchmark

        // Calculate the propagation rate
        // @see
        // https://github.com/google/benchmark/blob/main/docs/user_guide.md#custom-counters
        std::size_t total_tracks = 0u;
        for (auto _ : state) {
#pragma omp parallel for
            for (int i = 0; i < n_samples; ++i) {
                run_propagation((*tracks)[static_cast<std::size_t>(i)]);
            }
            total_tracks += static_cast<std::size_t>(n_samples);
        }

        // Report throughput
        state.counters["TracksPropagated"] = benchmark::Counter(
            static_cast<double>(total_tracks), benchmark::Counter::kIsRate);
    }
};

}  // namespace detray::benchmarks
