/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <vecmem/memory/binary_page_memory_resource.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "benchmark_propagator_cuda_kernel.hpp"
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

void fill_tracks(vecmem::vector<free_track_parameters> &tracks,
                 const unsigned int theta_steps, const unsigned int phi_steps) {
    // Set origin position of tracks
    const point3 ori{0., 0., 0.};

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.001 + itheta * (M_PI - 0.001) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            vector3 mom{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};
            mom = 10. * mom;

            // intialize a track
            free_track_parameters traj(ori, 0, mom, -1);

            tracks.push_back(traj);
        }
    }
}

static void BM_PROPAGATOR_CPU(benchmark::State &state) {

    // Create the toy geometry
    detector_host_type det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(host_mr, n_brl_layers,
                                                   n_edc_layers);

    // Set the magnetic field
    const vector3 B{0, 0, 2 * unit_constants::T};
    field_type B_field(B);

    // Create RK stepper
    rk_stepper_type s(B_field);

    // Create navigator
    navigator_host_type n(det);

    // Create propagator
    propagator_host_type p(std::move(s), std::move(n));

    for (auto _ : state) {

        state.PauseTiming();

        // Get tracks
        vecmem::vector<free_track_parameters> tracks(&host_mr);
        fill_tracks(tracks, state.range(0), state.range(0));

        state.ResumeTiming();

        for (auto &track : tracks) {

            // Create the propagator state
            propagator_host_type::state p_state(track);

            // Run propagation
            p.propagate(p_state);
        }
    }
}

static void BM_PROPAGATOR_CUDA(benchmark::State &state) {

    // Create the toy geometry
    detector_host_type det =
        create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                            vecmem::jagged_vector>(bp_mng_mr, n_brl_layers,
                                                   n_edc_layers);

    // Get detector data
    auto det_data = get_data(det);

    // vecmem copy helper object
    vecmem::cuda::copy copy;

    for (auto _ : state) {

        state.PauseTiming();

        // Get tracks
        vecmem::vector<free_track_parameters> tracks(&bp_mng_mr);
        fill_tracks(tracks, state.range(0), state.range(0));

        state.ResumeTiming();

        // Get tracks data
        auto tracks_data = vecmem::get_data(tracks);

        // Create navigator candidates buffer
        auto candidates_buffer =
            create_candidates_buffer(det, tracks.size(), dev_mr);
        copy.setup(candidates_buffer);

        // Run the propagator test for GPU device
        propagator_benchmark(det_data, tracks_data, candidates_buffer);
    }
}

BENCHMARK(BM_PROPAGATOR_CPU)->RangeMultiplier(2)->Range(8, 256);
BENCHMARK(BM_PROPAGATOR_CUDA)->RangeMultiplier(2)->Range(8, 256);

BENCHMARK_MAIN();