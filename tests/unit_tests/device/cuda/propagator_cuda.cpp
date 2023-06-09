/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "propagator_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

using namespace detray;

class CudaPropagatorConstBField : public ::testing::TestWithParam<vector3_t> {};

TEST_P(CudaPropagatorConstBField, propagator) {

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;

    // Set the magnetic field
    const vector3_t B = GetParam();

    // Create the toy geometry
    toy_cfg.bfield_vec(B);
    auto det = create_toy_geometry<const_bfield_bknd_t, host_container_types>(
        mng_mr, toy_cfg);

    // Create the vector of initial track parameterizations
    auto tracks_host = generate_tracks(&mng_mr);
    vecmem::vector<track_t> tracks_device(tracks_host, &mng_mr);

    // Host propagation
    auto&& [host_path_lengths, host_positions, host_jac_transports] =
        run_propagation_host(&mng_mr, det, tracks_host);

    // Device propagation
    auto&& [device_path_lengths, device_positions, device_jac_transports] =
        run_propagation_device<const_bfield_bknd_t>(&mng_mr, det, tracks_device,
                                                    host_positions);

    // Check the results
    compare_propagation_results(host_positions, device_positions,
                                host_path_lengths, device_path_lengths,
                                host_jac_transports, device_jac_transports);
}

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation1, CudaPropagatorConstBField,
                         ::testing::Values(vector3_t{0. * unit<scalar>::T,
                                                     0. * unit<scalar>::T,
                                                     2. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation2, CudaPropagatorConstBField,
                         ::testing::Values(vector3_t{0. * unit<scalar>::T,
                                                     1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation3, CudaPropagatorConstBField,
                         ::testing::Values(vector3_t{1. * unit<scalar>::T,
                                                     0. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation4, CudaPropagatorConstBField,
                         ::testing::Values(vector3_t{1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

/// This tests the device propagation in an inhomogenepus magnetic field
TEST(CudaPropagatorValidation5, inhomogeneous_bfield) {

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;

    // Create the toy geometry with inhomogeneous bfield from file
    auto det = create_toy_geometry<inhom_bfield_bknd_t, host_container_types>(
        mng_mr, toy_cfg);

    // Create the vector of initial track parameterizations
    auto tracks_host = generate_tracks(&mng_mr);
    vecmem::vector<track_t> tracks_device(tracks_host, &mng_mr);

    // Host propagation
    auto&& [host_path_lengths, host_positions, host_jac_transports] =
        run_propagation_host(&mng_mr, det, tracks_host);

    // Device propagation
    auto&& [device_path_lengths, device_positions, device_jac_transports] =
        run_propagation_device<inhom_bfield_bknd_t>(&mng_mr, det, tracks_device,
                                                    host_positions);

    // Check the results
    compare_propagation_results(host_positions, device_positions,
                                host_path_lengths, device_path_lengths,
                                host_jac_transports, device_jac_transports);
}
