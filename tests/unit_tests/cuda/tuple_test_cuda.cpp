/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <iostream>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "tuple_test_cuda_kernel.hpp"

TEST(tuple_test_cuda, tuple_test) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // mask_store-like tuple container
    vec_tuple<vecmem::vector, int, float, double> input_host(mng_mr);
    std::get<0>(input_host._tuple) = vecmem::vector<int>({1, 2, 3});
    std::get<1>(input_host._tuple) = vecmem::vector<float>({1.1, 5, 6, 3.4});
    std::get<2>(input_host._tuple) = vecmem::vector<double>({2.1, 2.2, 0.});

    vec_tuple_data input_data(input_host);

    // test std::make_tuple
    // std::tuple<int,int,int> test = std::make_tuple(1,1,1);

    tuple_test(input_data);
}