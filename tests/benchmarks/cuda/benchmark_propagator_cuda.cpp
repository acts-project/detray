/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "benchmark_propagator_cuda_kernel.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"

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
toy_det_config toy_cfg{4u, 7u};

void fill_tracks(vecmem::vector<free_track_parameters<transform3>> &tracks,
                 const std::size_t theta_steps, const std::size_t phi_steps) {
    // Set momentum of tracks
    const scalar mom_mag{10.f * unit<scalar>::GeV};

    // Iterate through uniformly distributed momentum directions
    for (auto traj : uniform_track_generator<free_track_parameters<transform3>>(
             phi_steps, theta_steps, mom_mag)) {
        tracks.push_back(traj);
    }
}

template <propagate_option opt>
static void BM_PROPAGATOR_CPU(benchmark::State &state) {

    // Create the toy geometry and bfield
    auto [det, names] = create_toy_geometry(host_mr, toy_cfg);
    vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
    auto bfield = bfield::create_const_field(B);

    // Create propagator
    propagation::config<scalar> cfg{};
    cfg.navigation.search_window = {3u, 3u};
    propagator_host_type p{cfg};

    std::size_t total_tracks = 0;

    for (auto _ : state) {

        // TODO: use fixture to build tracks
        state.PauseTiming();

        // Get tracks
        vecmem::vector<free_track_parameters<transform3>> tracks(&host_mr);
        fill_tracks(tracks, static_cast<std::size_t>(state.range(0)),
                    static_cast<std::size_t>(state.range(0)));

        total_tracks += tracks.size();

        state.ResumeTiming();

#pragma omp parallel for
        for (auto &track : tracks) {

            parameter_transporter<transform3>::state transporter_state{};
            pointwise_material_interactor<transform3>::state interactor_state{};
            parameter_resetter<transform3>::state resetter_state{};

            auto actor_states =
                tie(transporter_state, interactor_state, resetter_state);

            // Create the propagator state
            propagator_host_type::state p_state(track, bfield, det);

            // Run propagation
            if constexpr (opt == propagate_option::e_unsync) {
                p.propagate(p_state, actor_states);
            } else if constexpr (opt == propagate_option::e_sync) {
                p.propagate_sync(p_state, actor_states);
            }
        }
    }

    state.counters["TracksPropagated"] = benchmark::Counter(
        static_cast<double>(total_tracks), benchmark::Counter::kIsRate);
}

template <propagate_option opt>
static void BM_PROPAGATOR_CUDA(benchmark::State &state) {

    // Create the toy geometry
    auto [det, names] = create_toy_geometry(bp_mng_mr, toy_cfg);
    vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
    auto bfield = bfield::create_const_field(B);

    // Get detector data
    auto det_data = detray::get_data(det);

    // vecmem copy helper object
    vecmem::cuda::copy copy;

    std::size_t total_tracks = 0;

    for (auto _ : state) {

        state.PauseTiming();

        // Get tracks
        vecmem::vector<free_track_parameters<transform3>> tracks(&bp_mng_mr);
        fill_tracks(tracks, static_cast<std::size_t>(state.range(0)),
                    static_cast<std::size_t>(state.range(0)));

        total_tracks += tracks.size();

        state.ResumeTiming();

        // Get tracks data
        auto tracks_data = vecmem::get_data(tracks);

        // Create navigator candidates buffer
        auto candidates_buffer =
            create_candidates_buffer(det, tracks.size(), dev_mr, &mng_mr);
        copy.setup(candidates_buffer);

        // Run the propagator test for GPU device
        propagator_benchmark(det_data, bfield, tracks_data, candidates_buffer,
                             opt);
    }

    state.counters["TracksPropagated"] = benchmark::Counter(
        static_cast<double>(total_tracks), benchmark::Counter::kIsRate);
}

BENCHMARK_TEMPLATE(BM_PROPAGATOR_CPU, propagate_option::e_unsync)
    ->Name("CPU unsync propagation")
    ->RangeMultiplier(2)
    ->Range(8, 256);
BENCHMARK_TEMPLATE(BM_PROPAGATOR_CPU, propagate_option::e_sync)
    ->Name("CPU sync propagation")
    ->RangeMultiplier(2)
    ->Range(8, 256);
BENCHMARK_TEMPLATE(BM_PROPAGATOR_CUDA, propagate_option::e_unsync)
    ->Name("CUDA unsync propagation")
    ->RangeMultiplier(2)
    ->Range(8, 256);
BENCHMARK_TEMPLATE(BM_PROPAGATOR_CUDA, propagate_option::e_sync)
    ->Name("CUDA sync propagation")
    ->RangeMultiplier(2)
    ->Range(8, 256);

BENCHMARK_MAIN();
