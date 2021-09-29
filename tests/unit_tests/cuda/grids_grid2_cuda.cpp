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

#include "grids_grid2_cuda_kernel.hpp"

using namespace detray;

TEST(grids_cuda, grid2_replace_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // axis
    axis::regular<> xaxis{4, -1., 3.};
    axis::regular<> yaxis{6, 0., 6.};

    auto x_interval = (xaxis.max - xaxis.min) / xaxis.n_bins;
    auto y_interval = (yaxis.max - yaxis.min) / yaxis.n_bins;

    // declare host grid
    grid2<host_replace, axis::regular<>, axis::regular<>, serializer2> g2(
        std::move(xaxis), std::move(yaxis), mng_mr, test::point3{0, 0, 0});

    // pre-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {

            auto& data = g2.bin(i_x, i_y);

            EXPECT_EQ(data, g2.populator().kInvalid);
        }
    }

    // get grid_data
    grid2_data<device_replace, axis::regular<>, axis::regular<>, serializer2>
        g2_data(g2, mng_mr);

    // fill the grids
    grid_replace_test(g2_data);

    // post-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
            auto bin_id = i_x + i_y * xaxis.bins();
            auto& data = g2.bin(i_x, i_y);

            test::point3 tp({xaxis.min + bin_id * x_interval,
                             yaxis.min + bin_id * y_interval, 0.5});

            EXPECT_EQ(data, tp);
        }
    }
}

TEST(grids_cuda, grid2_complete_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // axis
    axis::regular<> xaxis{7, -1., 6.};
    axis::regular<> yaxis{3, 0., 3.};

    // declare grid
    grid2<host_complete, axis::regular<>, axis::regular<>, serializer2> g2(
        std::move(xaxis), std::move(yaxis), mng_mr, test::point3{0, 0, 0});

    // pre-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {

            auto& data = g2.bin(i_x, i_y);

            for (auto pt : data) {
                EXPECT_EQ(pt, g2.populator().kInvalid);
            }
        }
    }

    // get grid_data
    grid2_data<device_complete, axis::regular<>, axis::regular<>, serializer2>
        g2_data(g2, mng_mr);

    // fill the grid
    grid_complete_test(g2_data);

    auto x_interval = (xaxis.max - xaxis.min) / xaxis.n_bins;
    auto y_interval = (yaxis.max - yaxis.min) / yaxis.n_bins;

    // post-check
    for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
        for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {

            auto& data = g2.bin(i_x, i_y);

            for (int i_p = 0; i_p < data.size(); i_p++) {
                auto& pt = data[i_p];

                auto bin_id = i_x + i_y * xaxis.bins();
                auto gid = i_p + bin_id * data.size();

                test::point3 tp({xaxis.min + gid * x_interval,
                                 yaxis.min + gid * y_interval, 0.5});
                EXPECT_EQ(pt, tp);
            }
        }
    }
}

/**
 * This test demonstrates how to call grid buffer without calling host grid
 *object It is especially useful when you don't need to save the objects in host
 *side (e.g. internal spacepoint creation in traccc)
 **/
TEST(grids_cuda, grid2_buffer_attach_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    axis::regular<> xaxis{2, -1., 3.};
    axis::regular<> yaxis{2, 0., 6.};

    grid2_buffer<device_attach_populator<false, test::point3>, axis::regular<>,
                 axis::regular<>, serializer2>
        g2_buffer(xaxis, yaxis, {2, 5, 8, 10}, {100, 200, 300, 400}, mng_mr);

    // Check if the initialization work well
    // Non-zero starting size not working yet so initial argument for sizes is
    // ignored (acts-projects/vecmem#95)
    auto& ptr = g2_buffer._buffer.m_ptr;
    EXPECT_EQ(ptr[0].size(), 0);
    EXPECT_EQ(ptr[1].size(), 0);
    EXPECT_EQ(ptr[2].size(), 0);
    EXPECT_EQ(ptr[3].size(), 0);
    EXPECT_EQ(ptr[0].capacity(), 100);
    EXPECT_EQ(ptr[1].capacity(), 200);
    EXPECT_EQ(ptr[2].capacity(), 300);
    EXPECT_EQ(ptr[3].capacity(), 400);

    grid2_data<device_attach_populator<false, test::point3>, axis::regular<>,
               axis::regular<>, serializer2>
        g2_data(g2_buffer);

    // fill each bin with 10 points
    grid_attach_test(g2_data);

    // Check if each bin has 10 points
    EXPECT_EQ(ptr[0].size(), 100);
    EXPECT_EQ(ptr[1].size(), 100);
    EXPECT_EQ(ptr[2].size(), 100);
    EXPECT_EQ(ptr[3].size(), 100);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
