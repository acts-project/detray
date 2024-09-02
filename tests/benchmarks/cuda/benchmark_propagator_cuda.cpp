/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "benchmark_propagator_cuda_kernel.hpp"
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/common/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/binary_page_memory_resource.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Google include(s).
#include <benchmark/benchmark.h>

using namespace detray;

// VecMem memory resource(s)
vecmem::host_memory_resource host_mr;
vecmem::cuda::managed_memory_resource mng_mr;
vecmem::cuda::device_memory_resource dev_mr;
vecmem::binary_page_memory_resource bp_mng_mr(mng_mr);

// detector configuration
auto toy_cfg =
    toy_det_config{}.n_brl_layers(4u).n_edc_layers(7u).do_check(false);

void fill_tracks(vecmem::vector<free_track_parameters<algebra_t>> &tracks,
                 const std::size_t theta_steps, const std::size_t phi_steps) {
    // Set momentum of tracks
    const scalar mom_mag{10.f * unit<scalar>::GeV};

    // Iterate through uniformly distributed momentum directions
    for (auto traj : uniform_track_generator<free_track_parameters<algebra_t>>(
             phi_steps, theta_steps, mom_mag)) {
        tracks.push_back(traj);
    }
}

template <propagate_option opt>
static void BM_PROPAGATOR_CUDA(benchmark::State &state) {

    // Create the toy geometry
    auto [det, names] = build_toy_detector(bp_mng_mr, toy_cfg);
    test::vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
    auto bfield = bfield::create_const_field(B);

    // Get detector data
    auto det_data = detray::get_data(det);

    // vecmem copy helper object
    vecmem::cuda::copy copy;

    std::size_t total_tracks = 0;

    for (auto _ : state) {

        state.PauseTiming();

        // Get tracks
        vecmem::vector<free_track_parameters<algebra_t>> tracks(&bp_mng_mr);
        fill_tracks(tracks, static_cast<std::size_t>(state.range(0)),
                    static_cast<std::size_t>(state.range(0)));

        total_tracks += tracks.size();

        state.ResumeTiming();

        // Get tracks data
        auto tracks_data = vecmem::get_data(tracks);

        // Run the propagator test for GPU device
        propagator_benchmark(det_data, bfield, tracks_data, opt);
    }

    state.counters["TracksPropagated"] = benchmark::Counter(
        static_cast<double>(total_tracks), benchmark::Counter::kIsRate);
}

BENCHMARK_TEMPLATE(BM_PROPAGATOR_CUDA, propagate_option::e_unsync)
    ->Name("CUDA unsync propagation")
    ->RangeMultiplier(2)
    ->Range(8, 256);
BENCHMARK_TEMPLATE(BM_PROPAGATOR_CUDA, propagate_option::e_sync)
    ->Name("CUDA sync propagation")
    ->RangeMultiplier(2)
    ->Range(8, 256);

BENCHMARK_MAIN();
