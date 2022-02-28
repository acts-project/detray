/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "benchmark_propagator_cuda_kernel.hpp"
#include "vecmem/utils/cuda/copy.hpp"

using namespace detray;

// Helper object for performing memory copies.
vecmem::cuda::copy copy;

// VecMem memory resource(s)
vecmem::cuda::managed_memory_resource mng_mr;
vecmem::cuda::device_memory_resource dev_mr;

// Create the toy geometry
detector_host_type det =
    create_toy_geometry<darray, thrust::tuple, vecmem::vector,
                        vecmem::jagged_vector>(mng_mr, n_brl_layers,
                                               n_edc_layers);

// Create the vector of initial track parameters
vecmem::vector<free_track_parameters> tracks_host(&mng_mr);
vecmem::vector<free_track_parameters> tracks_device(&mng_mr);

// Set origin position of tracks
const point3 ori{0., 0., 0.};

void fill_tracks(vecmem::vector<free_track_parameters> &tracks) {
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
    for (auto _ : state) {
    }
}

static void BM_PROPAGATOR_CUDA(benchmark::State &state) {
    for (auto _ : state) {
    }
}

BENCHMARK(BM_PROPAGATOR_CPU);
BENCHMARK(BM_PROPAGATOR_CUDA);

BENCHMARK_MAIN();