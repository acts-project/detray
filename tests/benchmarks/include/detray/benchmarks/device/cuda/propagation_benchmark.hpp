/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

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
#include <random>
#include <string>

namespace detray {

// Define propagator type
template <typename algebra_t>
using default_chain = actor_chain<dtuple, parameter_transporter<algebra_t>,
                                  pointwise_material_interactor<algebra_t>,
                                  parameter_resetter<algebra_t>>;

template <typename metadata_t, typename bfield_t>
using cuda_propagator_type =
    propagator<rk_stepper<covfie::field_view<bfield_t>,
                          typename detector<metadata_t>::algebra_type>,
               navigator<detector<metadata_t>>,
               default_chain<typename detector<metadata_t>::algebra_type>>;

/// Launch the propagation kernelfor benchmarking
///
/// @param cfg the propagation configuration
/// @param det_view the detector vecmem view
/// @param field_data the magentic field view (maybe an empty field)
/// @param tracks_data the track collection view
/// @param navigation_cache_view the navigation cache vecemem view
/// @param opt which propagation to run (sync vs. unsync)
template <typename propagator_t>
void run_propagation_kernel(
    const propagation::config &cfg,
    typename propagator_t::detector_type::view_type det_view,
    typename propagator_t::stepper_type::magnetic_field_type field_data,
    vecmem::data::vector_view<
        free_track_parameters<typename propagator_t::algebra_type>>
        tracks_data,
    const int n_samples, const propagate_option opt);

/// Device Propagation becnhmark
template <typename propagator_t, typename bfield_bknd_t,
          propagate_option opt = propagate_option::e_sync>
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
                           const bfield_bknd_t *bfield) const {

        // Helper object for performing memory copies (to CUDA devices)
        vecmem::cuda::copy cuda_cpy;

        const int n_samples{m_cfg.benchmark().n_samples()};

        // Copy the track collection to device
        auto track_buffer =
            detray::get_buffer(vecmem::get_data(*tracks), *dev_mr, cuda_cpy);

        // Copy the detector to device and get its view
        auto det_buffer = detray::get_buffer(*det, *dev_mr, cuda_cpy);
        auto det_view = detray::get_data(det_buffer);

        std::size_t total_tracks = 0;
        for (auto _ : state) {
            // Launch the propagator test for GPU device
            run_propagation_kernel<propagator_t>(m_cfg.propagation(), det_view,
                                                 *bfield, track_buffer,
                                                 n_samples, opt);

            total_tracks += n_samples;
        }
        // Report throughput
        state.counters["TracksPropagated"] = benchmark::Counter(
            static_cast<double>(total_tracks), benchmark::Counter::kIsRate);
    }
};

}  // namespace detray
