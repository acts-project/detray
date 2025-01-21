/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s).
#include "detray/test/utils/types.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/benchmarks/propagation_benchmark_config.hpp"
#include "detray/benchmarks/propagation_benchmark_utils.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <string>

namespace detray::benchmarks {

// Define propagator type
template <concepts::algebra algebra_t>
using empty_chain = actor_chain<>;

template <concepts::algebra algebra_t>
using default_chain = actor_chain<parameter_transporter<algebra_t>,
                                  pointwise_material_interactor<algebra_t>,
                                  parameter_resetter<algebra_t>>;

using const_field_t = bfield::const_bknd_t<test::scalar>;

template <typename metadata_t, typename bfield_t,
          template <typename> class actor_chain_t>
using cuda_propagator_type =
    propagator<rk_stepper<covfie::field_view<bfield_t>,
                          typename detector<metadata_t>::algebra_type>,
               navigator<detector<metadata_t>>,
               actor_chain_t<typename detector<metadata_t>::algebra_type>>;

/// Launch the propagation kernelfor benchmarking
///
/// @param cfg the propagation configuration
/// @param det_view the detector vecmem view
/// @param field_data the magentic field view (maybe an empty field)
/// @param tracks_data the track collection view
/// @param navigation_cache_view the navigation cache vecemem view
/// @param opt which propagation to run (sync vs. unsync)
template <typename propagator_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
void run_propagation_kernel(
    const propagation::config &,
    typename propagator_t::detector_type::view_type,
    typename propagator_t::stepper_type::magnetic_field_type,
    typename propagator_t::actor_chain_type::state_tuple *,
    vecmem::data::vector_view<
        free_track_parameters<typename propagator_t::algebra_type>>,
    const int);

/// Allocate actor state blueprint on device
/// @note This only works if each actor state in the tuple is essentially POD
template <typename propagator_t>
typename propagator_t::actor_chain_type::state_tuple *setup_actor_states(
    typename propagator_t::actor_chain_type::state_tuple *);

/// Release actor state blueprint
template <typename propagator_t>
void release_actor_states(
    typename propagator_t::actor_chain_type::state_tuple *);

/// Device Propagation becnhmark
template <typename propagator_t, typename bfield_bknd_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
struct cuda_propagation_bm : public benchmark_base {
    /// Detector dependent types
    using algebra_t = typename propagator_t::detector_type::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    /// Local configuration type
    using configuration = propagation_benchmark_config;

    /// The benchmark configuration
    configuration m_cfg{};

    /// Default construction
    cuda_propagation_bm() = default;

    /// Construct from an externally provided configuration @param cfg
    explicit cuda_propagation_bm(const configuration &cfg) : m_cfg{cfg} {}

    /// @return the benchmark configuration
    configuration &config() { return m_cfg; }

    /// Prepare data and run benchmark loop
    inline void operator()(::benchmark::State &state,
                           vecmem::memory_resource *dev_mr,
                           dvector<free_track_parameters<algebra_t>> *tracks,
                           const typename propagator_t::detector_type *det,
                           const bfield_bknd_t *bfield,
                           typename propagator_t::actor_chain_type::state_tuple
                               *input_actor_states) const {

        assert(dev_mr != nullptr);
        assert(tracks != nullptr);
        assert(det != nullptr);
        assert(bfield != nullptr);
        assert(input_actor_states != nullptr);

        // Helper object for performing memory copies (to CUDA devices)
        vecmem::cuda::copy cuda_cpy;

        const int n_samples{m_cfg.benchmark().n_samples()};
        const int n_warmup{m_cfg.benchmark().n_warmup()};

        assert(static_cast<std::size_t>(n_samples) <= tracks->size());

        // Copy the track collection to device
        auto track_buffer =
            detray::get_buffer(vecmem::get_data(*tracks), *dev_mr, cuda_cpy);

        // Copy the detector to device and get its view
        auto det_buffer = detray::get_buffer(*det, *dev_mr, cuda_cpy);
        auto det_view = detray::get_data(det_buffer);

        // Copy blueprint actor states to device
        auto *device_actor_state_ptr =
            setup_actor_states<propagator_t>(input_actor_states);

        // Do a small warm up run
        if (m_cfg.benchmark().do_warmup()) {
            auto warmup_track_buffer = detray::get_buffer(
                vecmem::get_data(*tracks), *dev_mr, cuda_cpy);

            run_propagation_kernel<propagator_t, kOPT>(
                m_cfg.propagation(), det_view, *bfield, device_actor_state_ptr,
                warmup_track_buffer, math::min(n_warmup, n_samples));
        } else {
            std::cout << "WARNING: Running CUDA benchmarks without warmup is "
                         "not recommended"
                      << std::endl;
        }

        // Calculate the propagation rate
        // @see
        // https://github.com/google/benchmark/blob/main/docs/user_guide.md#custom-counters
        std::size_t total_tracks = 0u;
        for (auto _ : state) {
            // Launch the propagator test for GPU device
            run_propagation_kernel<propagator_t, kOPT>(
                m_cfg.propagation(), det_view, *bfield, device_actor_state_ptr,
                track_buffer, n_samples);

            total_tracks += static_cast<std::size_t>(n_samples);
        }

        // Report throughput
        state.counters["TracksPropagated"] = benchmark::Counter(
            static_cast<double>(total_tracks), benchmark::Counter::kIsRate);

        release_actor_states<propagator_t>(device_actor_state_ptr);
    }
};

}  // namespace detray::benchmarks
