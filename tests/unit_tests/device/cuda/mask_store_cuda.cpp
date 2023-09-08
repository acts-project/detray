/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <cstdlib>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "mask_store_cuda_kernel.hpp"

TEST(mask_store_cuda, mask_store) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // Types must be sorted according to their id (here: masks/mask_identifier)
    host_store_type store(mng_mr);

    ASSERT_TRUE(store.template empty<e_annulus2>());
    ASSERT_TRUE(store.template empty<e_cylinder2>());
    ASSERT_TRUE(store.template empty<e_rectangle2>());
    ASSERT_TRUE(store.template empty<e_ring2>());
    ASSERT_TRUE(store.template empty<e_single3>());
    ASSERT_TRUE(store.template empty<e_trapezoid2>());

    store.template emplace_back<e_rectangle2>(empty_context{}, 0u, 1.0f, 2.0f);
    store.template emplace_back<e_trapezoid2>(empty_context{}, 0u, 0.5f, 1.5f,
                                              4.0f);
    store.template emplace_back<e_ring2>(empty_context{}, 0u, 1.0f, 10.0f);
    store.template emplace_back<e_cylinder2>(empty_context{}, 0u, 1.f, 0.5f,
                                             2.0f);
    store.template emplace_back<e_annulus2>(empty_context{}, 0u, 1.f, 2.f, 3.f,
                                            4.f, 5.f, 6.f, 7.f);

    ASSERT_FALSE(store.empty<e_annulus2>());
    ASSERT_FALSE(store.empty<e_cylinder2>());
    ASSERT_FALSE(store.empty<e_rectangle2>());
    ASSERT_FALSE(store.empty<e_ring2>());
    ASSERT_TRUE(store.empty<e_single3>());
    ASSERT_FALSE(store.empty<e_trapezoid2>());

    /** Generate random points for test **/
    vecmem::vector<point3> input_point3(n_points, &mng_mr);

    for (unsigned int i = 0u; i < n_points; i++) {
        point3 rand_point3 = {static_cast<scalar>(rand() % 100) / 10.f,
                              static_cast<scalar>(rand() % 100) / 10.f,
                              static_cast<scalar>(rand() % 100) / 10.f};
        input_point3[i] = rand_point3;
    }

    /** host output for intersection status **/
    vecmem::jagged_vector<bool> output_host(5, &mng_mr);

    /** get mask objects **/
    const auto& rectangle_mask = store.get<e_rectangle2>()[0];
    const auto& trapezoid_mask = store.get<e_trapezoid2>()[0];
    const auto& ring_mask = store.get<e_ring2>()[0];
    const auto& cylinder_mask = store.get<e_cylinder2>()[0];
    const auto& annulus_mask = store.get<e_annulus2>()[0];

    /** get host results from is_inside function **/
    for (unsigned int i = 0u; i < n_points; i++) {
        output_host[0].push_back(rectangle_mask.is_inside(input_point3[i]));
        output_host[1].push_back(trapezoid_mask.is_inside(input_point3[i]));
        output_host[2].push_back(ring_mask.is_inside(input_point3[i]));
        output_host[3].push_back(cylinder_mask.is_inside(input_point3[i]));
        output_host[4].push_back(annulus_mask.is_inside(input_point3[i]));
    }

    /** Helper object for performing memory copies. **/
    vecmem::cuda::copy copy;

    /** device output for intersection status **/
    vecmem::data::jagged_vector_buffer<bool> output_buffer(
        {n_points, n_points, n_points, n_points, n_points}, mng_mr, nullptr,
        vecmem::data::buffer_type::resizable);

    copy.setup(output_buffer);

    auto input_point3_data = vecmem::get_data(input_point3);
    // auto input_point3_data = vecmem::get_data(input_point3);
    auto store_data = get_data(store);

    /** run the kernel **/
    mask_test(store_data, input_point3_data, output_buffer);

    vecmem::jagged_vector<bool> output_device(&mng_mr);
    copy(output_buffer, output_device);

    /** Compare the values **/
    for (unsigned int i = 0u; i < n_points; i++) {
        ASSERT_EQ(output_host[0][i], output_device[0][i]);
        ASSERT_EQ(output_host[0][i], output_device[0][i]);
        ASSERT_EQ(output_host[2][i], output_device[2][i]);
        ASSERT_EQ(output_host[3][i], output_device[3][i]);
        ASSERT_EQ(output_host[4][i], output_device[4][i]);
    }
}
