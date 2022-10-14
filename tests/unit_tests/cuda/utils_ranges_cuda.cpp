/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "utils_ranges_cuda_kernel.hpp"
#include "vecmem/utils/cuda/copy.hpp"

using namespace detray;

// This tests the single value view
TEST(utils_ranges_cuda, single) {

    dindex value{251UL};
    dindex check{std::numeric_limits<dindex>::max()};

    // Run the test code
    test_single(value, check);

    // Check result value
    ASSERT_EQ(value, check);
}

// This tests the iota range generator
TEST(utils_ranges_cuda, iota) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input reference vector for test
    vecmem::vector<dindex> reference({2, 3, 4, 5, 6});

    // Const range for enumeration test
    const darray<dindex, 2> range = {2, 7};

    // Output vector buffer for enumeration test
    vecmem::data::vector_buffer<dindex> check_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(range[1] -
                                                                    range[0]),
        0, managed_resource);
    copy.setup(check_buffer);

    // Run test function
    test_iota(range, check_buffer);

    // Copy vector buffer to output vector
    vecmem::vector<dindex> check{&managed_resource};
    copy(check_buffer, check);

    // Check the result
    ASSERT_EQ(check, reference);
}

// This tests the convenience enumeration function
TEST(utils_ranges_cuda, enumerate) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input vector sequence for test
    vecmem::vector<uint_holder> seq({{0}, {1}, {2}, {3}, {4}, {5}},
                                    &managed_resource);

    // Get vector_data object
    auto seq_data = vecmem::get_data(seq);

    // Output vector buffer for enumeration test
    vecmem::data::vector_buffer<dindex> idx_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
        0, managed_resource);
    copy.setup(idx_buffer);

    vecmem::data::vector_buffer<dindex> value_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
        0, managed_resource);
    copy.setup(value_buffer);

    // Run test function
    test_enumerate(seq_data, idx_buffer, value_buffer);

    // Copy vector buffer to output vector
    vecmem::vector<dindex> idx_vec{&managed_resource};
    copy(idx_buffer, idx_vec);

    vecmem::vector<dindex> value_vec{&managed_resource};
    copy(value_buffer, value_vec);

    // Check the result
    for (std::size_t i = 0; i < idx_vec.size(); i++) {
        ASSERT_EQ(idx_vec[i], value_vec[i]);
    }
}

// This tests the convenience pick function
TEST(utils_ranges_cuda, pick) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input vector sequence for test
    vecmem::vector<uint_holder> seq({{0}, {1}, {2}, {3}, {4}, {5}},
                                    &managed_resource);
    // Input index sequence for test
    vecmem::vector<dindex> idx({0, 2, 4, 5}, &managed_resource);

    // Get vector_data object
    auto seq_data = vecmem::get_data(seq);
    auto idx_data = vecmem::get_data(idx);

    // Output vector buffer for enumeration test
    vecmem::data::vector_buffer<dindex> idx_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
        0, managed_resource);
    copy.setup(idx_buffer);

    vecmem::data::vector_buffer<dindex> value_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
        0, managed_resource);
    copy.setup(value_buffer);

    // Run test function
    test_pick(seq_data, idx_data, idx_buffer, value_buffer);

    // Copy vector buffer to output vector
    vecmem::vector<dindex> idx_vec{&managed_resource};
    copy(idx_buffer, idx_vec);

    vecmem::vector<dindex> value_vec{&managed_resource};
    copy(value_buffer, value_vec);

    // Check the result
    for (std::size_t i = 0; i < idx_vec.size(); i++) {
        ASSERT_EQ(idx[i], idx_vec[i]);
        ASSERT_EQ(seq[idx_vec[i]].ui, value_vec[i]);
    }
}

// This tests the convenience enumeration function
TEST(utils_ranges_cuda, join) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input vector sequence for test
    vecmem::vector<uint_holder> seq_1({{0}, {1}, {2}, {3}, {4}, {5}},
                                      &managed_resource);
    vecmem::vector<uint_holder> seq_2({{2}, {0}, {9}, {4}, {15}},
                                      &managed_resource);
    // Get vector_data object
    auto seq_data_1 = vecmem::get_data(seq_1);
    auto seq_data_2 = vecmem::get_data(seq_2);

    // Output vector buffer for join test
    vecmem::data::vector_buffer<dindex> value_buffer(
        static_cast<vecmem::data::vector_buffer<dindex>::size_type>(
            seq_1.size() + seq_2.size()),
        0, managed_resource);
    copy.setup(value_buffer);

    // Run test function
    test_join(seq_data_1, seq_data_2, value_buffer);

    // Copy vector buffer to output vector
    vecmem::vector<dindex> value_vec{&managed_resource};
    copy(value_buffer, value_vec);

    // First sequence
    for (std::size_t i = 0; i < seq_1.size(); i++) {
        ASSERT_EQ(seq_1[i].ui, value_vec[i]);
    }
    // Second sequence
    for (std::size_t i = 0; i < seq_1.size(); i++) {
        ASSERT_EQ(seq_2[i].ui, value_vec[i + seq_1.size()]);
    }
}

// This tests the subrange view
TEST(utils_ranges_cuda, subrange) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource managed_resource;

    // Input vector sequence for test
    vecmem::vector<int> seq({0, 1, 2, 3, 4, 5}, &managed_resource);
    // Get vector_data object
    auto seq_data = vecmem::get_data(seq);

    // Begin and end index for iteration
    const std::size_t begin = 1;
    const std::size_t end = 4;

    // Output vector buffer for iteration test
    vecmem::data::vector_buffer<int> check_buffer(
        static_cast<vecmem::data::vector_buffer<int>::size_type>(end - begin),
        0, managed_resource);
    copy.setup(check_buffer);

    // Run test function
    test_subrange(seq_data, check_buffer, begin, end);

    // Copy vector buffer to output vector
    vecmem::vector<int> check{&managed_resource};
    copy(check_buffer, check);

    // Check the result
    ASSERT_EQ(check[0], 1);
    ASSERT_EQ(check[1], 2);
    ASSERT_EQ(check[2], 3);
}
