/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "propagator_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include
#include <gtest/gtest.h>

using namespace detray;

class CudaPropConstBFieldMng : public ::testing::TestWithParam<vector3_t> {};

TEST_P(CudaPropConstBFieldMng, propagator) {

    // VecMem memory resource(s)
    vecmem::cuda::managed_memory_resource mng_mr;

    // Set the magnetic field
    const vector3_t B = GetParam();

    // Create the toy geometry
    toy_cfg.bfield_vec(B);
    auto det = create_toy_geometry<const_bfield_bknd_t>(mng_mr, toy_cfg);

    run_propagation_test(&mng_mr, det, detray::get_data(det));
}

class CudaPropConstBFieldCpy : public ::testing::TestWithParam<vector3_t> {};

TEST_P(CudaPropConstBFieldCpy, propagator) {

    // VecMem memory resource(s)
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    vecmem::cuda::copy cuda_cpy;

    // Set the magnetic field
    const vector3_t B = GetParam();

    // Create the toy geometry
    toy_cfg.bfield_vec(B);
    auto det = create_toy_geometry<const_bfield_bknd_t>(host_mr, toy_cfg);

    auto det_buff = detray::get_buffer<covfie::field<const_bfield_bknd_t>>(
        det, dev_mr, cuda_cpy);

    run_propagation_test(&mng_mr, det, detray::get_data(det_buff));
}

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation1, CudaPropConstBFieldMng,
                         ::testing::Values(vector3_t{0. * unit<scalar>::T,
                                                     0. * unit<scalar>::T,
                                                     2. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation2, CudaPropConstBFieldMng,
                         ::testing::Values(vector3_t{0. * unit<scalar>::T,
                                                     1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation3, CudaPropConstBFieldMng,
                         ::testing::Values(vector3_t{1. * unit<scalar>::T,
                                                     0. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation4, CudaPropConstBFieldMng,
                         ::testing::Values(vector3_t{1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation5, CudaPropConstBFieldCpy,
                         ::testing::Values(vector3_t{0. * unit<scalar>::T,
                                                     0. * unit<scalar>::T,
                                                     2. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation6, CudaPropConstBFieldCpy,
                         ::testing::Values(vector3_t{0. * unit<scalar>::T,
                                                     1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation7, CudaPropConstBFieldCpy,
                         ::testing::Values(vector3_t{1. * unit<scalar>::T,
                                                     0. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

INSTANTIATE_TEST_SUITE_P(CudaPropagatorValidation8, CudaPropConstBFieldCpy,
                         ::testing::Values(vector3_t{1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T,
                                                     1. * unit<scalar>::T}));

/// This tests the device propagation in an inhomogenepus magnetic field
TEST(CudaPropagatorValidation10, inhomogeneous_bfield_cpy) {

    // VecMem memory resource(s)
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    vecmem::cuda::copy cuda_cpy;

    // Create the toy geometry with inhomogeneous bfield from file
    auto det = create_toy_geometry<inhom_bfield_bknd_t>(host_mr, toy_cfg);

    auto det_buff = detray::get_buffer<covfie::field<inhom_bfield_bknd_cuda_t>>(
        det, dev_mr, cuda_cpy);

    run_propagation_test(&mng_mr, det, detray::get_data(det_buff));
}
