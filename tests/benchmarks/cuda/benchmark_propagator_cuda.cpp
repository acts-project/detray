/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <covfie/core/field.hpp>
#include <vecmem/memory/binary_page_memory_resource.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "benchmark_propagator_cuda_kernel.hpp"
#include "detray/simulation/track_generators.hpp"
#include "vecmem/utils/cuda/copy.hpp"

using namespace detray;

// VecMem memory resource(s)
vecmem::host_memory_resource host_mr;
vecmem::cuda::managed_memory_resource mng_mr;
vecmem::cuda::device_memory_resource dev_mr;
vecmem::binary_page_memory_resource bp_mng_mr(mng_mr);

// detector configuration
constexpr std::size_t n_brl_layers = 4;
constexpr std::size_t n_edc_layers = 7;

void fill_tracks(vecmem::vector<free_track_parameters<transform3>> &tracks,
                 const unsigned int theta_steps, const unsigned int phi_steps) {
    // Set origin position of tracks
    const point3 ori{0., 0., 0.};
    const scalar mom_mag = 10. * unit<scalar>::GeV;

    // Iterate through uniformly distributed momentum directions
    for (auto traj : uniform_track_generator<free_track_parameters<transform3>>(
             theta_steps, phi_steps, ori, mom_mag)) {
        tracks.push_back(traj);
    }
}

static void BM_PROPAGATOR_CPU(benchmark::State &state) {

    // Create the toy geometry
    detector_host_type det = create_toy_geometry<host_container_types>(
        host_mr,
        field_type(field_type::backend_t::configuration_t{
            0.f, 0.f, 2.f * unit<scalar>::T}),
        n_brl_layers, n_edc_layers);

    // Create RK stepper
    rk_stepper_type s;

    // Create navigator
    navigator_host_type n;

    // Create propagator
    propagator_host_type p(std::move(s), std::move(n));

    for (auto _ : state) {

        // TODO: use fixture to build tracks
        state.PauseTiming();

        // Get tracks
        vecmem::vector<free_track_parameters<transform3>> tracks(&host_mr);
        fill_tracks(tracks, state.range(0), state.range(0));

        state.ResumeTiming();

        for (auto &track : tracks) {

            parameter_transporter<transform3>::state transporter_state{};
            pointwise_material_interactor<transform3>::state interactor_state{};
            parameter_resetter<transform3>::state resetter_state{};

            auto actor_states = thrust::tie(transporter_state, interactor_state,
                                            resetter_state);

            // Create the propagator state
            propagator_host_type::state p_state(track, det.get_bfield(), det);

            // Run propagation
            p.propagate(p_state, actor_states);
        }
    }
}

static void BM_PROPAGATOR_CUDA(benchmark::State &state) {

    // Create the toy geometry
    detector_host_type det = create_toy_geometry<host_container_types>(
        bp_mng_mr, n_brl_layers, n_edc_layers);

    // Get detector data
    auto det_data = get_data(det);

    // vecmem copy helper object
    vecmem::cuda::copy copy;

    for (auto _ : state) {

        state.PauseTiming();

        // Get tracks
        vecmem::vector<free_track_parameters<transform3>> tracks(&bp_mng_mr);
        fill_tracks(tracks, state.range(0), state.range(0));

        state.ResumeTiming();

        // Get tracks data
        auto tracks_data = vecmem::get_data(tracks);

        // Create navigator candidates buffer
        auto candidates_buffer =
            create_candidates_buffer(det, tracks.size(), dev_mr, &mng_mr);
        copy.setup(candidates_buffer);

        // Run the propagator test for GPU device
        propagator_benchmark(det_data, tracks_data, candidates_buffer);
    }
}

BENCHMARK(BM_PROPAGATOR_CPU)->RangeMultiplier(2)->Range(8, 256);
BENCHMARK(BM_PROPAGATOR_CUDA)->RangeMultiplier(2)->Range(8, 256);

BENCHMARK_MAIN();
