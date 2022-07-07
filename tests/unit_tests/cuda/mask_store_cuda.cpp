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
    tuple_vector_container<thrust::tuple, dvector, mask_ids, rectangle,
                           trapezoid, ring, cylinder, single, annulus>
        store(mng_mr);

    ASSERT_TRUE(store.template empty<e_annulus2>());
    ASSERT_TRUE(store.template empty<e_cylinder3>());
    ASSERT_TRUE(store.template empty<e_rectangle2>());
    ASSERT_TRUE(store.template empty<e_ring2>());
    ASSERT_TRUE(store.template empty<e_single3>());
    ASSERT_TRUE(store.template empty<e_trapezoid2>());

    store.template add_value<e_rectangle2>(1.0, 2.0, 0);
    store.template add_value<e_trapezoid2>(0.5, 1.5, 4.0, 0);
    store.template add_value<e_ring2>(1.0, 10.0, 0);
    store.template add_value<e_cylinder3>(1., 0.5, 2.0, 0);
    // store.group<e_single3>().push_back(single{3.0,6.0});
    store.template add_value<e_annulus2>(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0);

    ASSERT_FALSE(store.empty<e_annulus2>());
    ASSERT_FALSE(store.empty<e_cylinder3>());
    ASSERT_FALSE(store.empty<e_rectangle2>());
    ASSERT_FALSE(store.empty<e_ring2>());
    // ASSERT_FALSE(store.template empty<e_single3>());
    ASSERT_FALSE(store.empty<e_trapezoid2>());

    /** Generate random points for test **/
    vecmem::vector<point2> input_point2(n_points, &mng_mr);
    vecmem::vector<point3> input_point3(n_points, &mng_mr);

    for (int i = 0; i < n_points; i++) {
        point2 rand_point2 = {static_cast<scalar>(rand() % 100 / 10.),
                              static_cast<scalar>(rand() % 100 / 10.)};
        point3 rand_point3 = {static_cast<scalar>(rand() % 100 / 10.),
                              static_cast<scalar>(rand() % 100 / 10.),
                              static_cast<scalar>(rand() % 100 / 10.)};

        input_point2[i] = rand_point2;
        input_point3[i] = rand_point3;
    }

    /** host output for intersection status **/
    vecmem::jagged_vector<intersection::status> output_host(5, &mng_mr);

    /** get mask objects **/
    const auto& rectangle_mask = store.group<e_rectangle2>()[0];
    const auto& trapezoid_mask = store.group<e_trapezoid2>()[0];
    const auto& ring_mask = store.group<e_ring2>()[0];
    const auto& cylinder_mask = store.group<e_cylinder3>()[0];
    const auto& annulus_mask = store.group<e_annulus2>()[0];

    /** get host results from is_inside function **/
    for (int i = 0; i < n_points; i++) {
        output_host[0].push_back(rectangle_mask.is_inside(input_point2[i]));
        output_host[1].push_back(trapezoid_mask.is_inside(input_point2[i]));
        output_host[2].push_back(ring_mask.is_inside(input_point2[i]));
        output_host[3].push_back(cylinder_mask.is_inside(input_point3[i]));
        output_host[4].push_back(annulus_mask.is_inside(input_point2[i]));
    }

    /** Helper object for performing memory copies. **/
    vecmem::cuda::copy copy;

    /** device output for intersection status **/
    vecmem::data::jagged_vector_buffer<intersection::status> output_buffer(
        {0, 0, 0, 0, 0}, {n_points, n_points, n_points, n_points, n_points},
        mng_mr);

    copy.setup(output_buffer);

    auto input_point2_data = vecmem::get_data(input_point2);
    auto input_point3_data = vecmem::get_data(input_point3);
    auto store_data = get_data(store);

    /** run the kernel **/
    mask_test(store_data, input_point2_data, input_point3_data, output_buffer);

    vecmem::jagged_vector<intersection::status> output_device(&mng_mr);
    copy(output_buffer, output_device);

    /** Compare the values **/
    for (int i = 0; i < n_points; i++) {
        ASSERT_EQ(output_host[0][i], output_device[0][i]);
        ASSERT_EQ(output_host[0][i], output_device[0][i]);
        ASSERT_EQ(output_host[2][i], output_device[2][i]);
        ASSERT_EQ(output_host[3][i], output_device[3][i]);
        ASSERT_EQ(output_host[4][i], output_device[4][i]);
    }
}
