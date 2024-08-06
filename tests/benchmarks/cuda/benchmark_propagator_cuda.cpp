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
auto toy_cfg = toy_det_config{}
                   .n_brl_layers(4u)
                   .n_edc_layers(7u)
                   .do_check(false)
                   .use_material_maps(false);

void fill_tracks(vecmem::vector<free_track_parameters<algebra_t>> &tracks,
                 const std::size_t theta_steps, const std::size_t phi_steps) {
    using scalar_t = dscalar<algebra_t>;
    using uniform_gen_t =
        detail::random_numbers<scalar_t,
                               std::uniform_real_distribution<scalar_t>>;
    using track_generator_t =
        random_track_generator<free_track_parameters<algebra_t>, uniform_gen_t>;

    track_generator_t::configuration trk_cfg{};
    trk_cfg.n_tracks(theta_steps * phi_steps)
        .randomize_charge(true)
        .eta_range(-4.f, 4.f)
        .pT_range(1.f * unit<scalar_t>::GeV, 100.f * unit<scalar_t>::GeV);

    // Iterate through uniformly distributed momentum directions
    for (auto traj : track_generator_t(trk_cfg)) {
        tracks.push_back(traj);
    }
}

template <propagate_option opt>
static void BM_PROPAGATOR_CPU(benchmark::State &state) {

    // Create the toy geometry and bfield
    auto [det, names] = build_toy_detector(host_mr, toy_cfg);
    test::vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
    auto bfield = bfield::create_const_field(B);

    // Create propagator
    propagation::config cfg{};
    cfg.navigation.search_window = {3u, 3u};
    propagator_host_type p{cfg};

    std::size_t total_tracks = 0;

    // Get tracks
    vecmem::vector<free_track_parameters<algebra_t>> tracks(&host_mr);
    fill_tracks(tracks, static_cast<std::size_t>(state.range(0)),
                static_cast<std::size_t>(state.range(0)));

    total_tracks += tracks.size();

    for (auto _ : state) {

#pragma omp parallel for
        for (auto &track : tracks) {

            parameter_transporter<algebra_t>::state transporter_state{};
            pointwise_material_interactor<algebra_t>::state interactor_state{};
            parameter_resetter<algebra_t>::state resetter_state{};

            auto actor_states =
                tie(transporter_state, interactor_state, resetter_state);

            // Create the propagator state
            propagator_host_type::state p_state(track, bfield, det);

            // Run propagation
            if constexpr (opt == propagate_option::e_unsync) {
                ::benchmark::DoNotOptimize(p.propagate(p_state, actor_states));
            } else if constexpr (opt == propagate_option::e_sync) {
                ::benchmark::DoNotOptimize(
                    p.propagate_sync(p_state, actor_states));
            }
        }
    }

    state.counters["TracksPropagated"] = benchmark::Counter(
        static_cast<double>(total_tracks), benchmark::Counter::kIsRate);
}

template <propagate_option opt>
static void BM_PROPAGATOR_CUDA(benchmark::State &state) {

    // Create the toy geometry
    auto [det, names] = build_toy_detector(bp_mng_mr, toy_cfg);
    test::vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
    auto bfield = bfield::create_const_field(B);

    // vecmem copy helper object
    vecmem::cuda::copy copy;

    // Copy the detector to device and get its view
    auto det_buffer = detray::get_buffer(det, dev_mr, copy);
    auto det_view = detray::get_data(det_buffer);

    std::size_t total_tracks = 0;

    // Get tracks
    vecmem::vector<free_track_parameters<algebra_t>> tracks(&bp_mng_mr);
    fill_tracks(tracks, static_cast<std::size_t>(state.range(0)),
                static_cast<std::size_t>(state.range(0)));

    total_tracks += tracks.size();

    // Get tracks data
    auto track_buffer =
        detray::get_buffer(vecmem::get_data(tracks), dev_mr, copy);

    // Create navigator candidates buffer
    auto candidates_buffer =
        create_candidates_buffer(det, tracks.size(), dev_mr, &mng_mr);
    copy.setup(candidates_buffer);

    for (auto _ : state) {
        // Run the propagator test for GPU device
        propagator_benchmark(det_view, bfield, track_buffer, candidates_buffer,
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
