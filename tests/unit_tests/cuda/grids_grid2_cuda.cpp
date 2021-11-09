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
#include "vecmem/utils/cuda/copy.hpp"

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
    host_grid2_replace g2(std::move(xaxis), std::move(yaxis), mng_mr,
                          test::point3{0, 0, 0});

    // pre-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {

            const auto& data = g2.bin(i_x, i_y);

            EXPECT_EQ(data, g2.populator().m_invalid);
        }
    }

    // get grid_data
    auto g2_data = get_data(g2, mng_mr);

    // fill the grids
    grid_replace_test(g2_data);

    // post-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
            auto bin_id = i_x + i_y * xaxis.bins();
            const auto& data = g2.bin(i_x, i_y);

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
    host_grid2_complete g2(std::move(xaxis), std::move(yaxis), mng_mr,
                           test::point3{0, 0, 0});

    // pre-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {

            const auto& data = g2.bin(i_x, i_y);

            for (auto pt : data) {
                EXPECT_EQ(pt, g2.populator().m_invalid);
            }
        }
    }

    // get grid_data
    auto g2_data = get_data(g2, mng_mr);

    // fill the grid
    grid_complete_test(g2_data);

    auto x_interval = (xaxis.max - xaxis.min) / xaxis.n_bins;
    auto y_interval = (yaxis.max - yaxis.min) / yaxis.n_bins;

    // post-check
    for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
        for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {

            const auto& data = g2.bin(i_x, i_y);

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

TEST(grids_cuda, grid2_attach_populator) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    axis::circular<> xaxis{65, -M_PI, M_PI};
    axis::regular<> yaxis{2, 0., 6.};

    auto x_interval = (xaxis.max - xaxis.min) / xaxis.n_bins;
    auto y_interval = (yaxis.max - yaxis.min) / yaxis.n_bins;

    host_grid2_attach g2(xaxis, yaxis, mng_mr, test::point3{0, 0, 0});

    for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
        for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {

            for (int i_p = 0; i_p < 100; i_p++) {

                auto bin_id = i_x + i_y * xaxis.bins();
                auto gid = i_p + bin_id * 100;

                test::point3 tp({xaxis.min + gid * x_interval,
                                 yaxis.min + gid * y_interval, 0.5});
                g2.populate(i_x, i_y, std::move(tp));
            }
        }
    }

    // get grid_data
    auto g2_data = get_data(g2, mng_mr);

    // Read the grid
    grid_attach_read_test(g2_data);
}

/// This test demonstrates how to call grid buffer without calling host grid
/// object It is especially useful when you don't need to save the objects in
/// host side (e.g. internal spacepoint creation in traccc)
TEST(grids_cuda, grid2_buffer_attach_populator) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    axis::circular<> xaxis{2, -1., 3.};
    axis::regular<> yaxis{2, 0., 6.};

    grid2_buffer<host_grid2_attach> g2_buffer(xaxis, yaxis, {2, 5, 8, 10},
                                              {100, 200, 300, 400}, mng_mr);
    copy.setup(g2_buffer._buffer);

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

    // fill each bin with 100 points
    grid_attach_fill_test(g2_buffer);

    host_grid2_attach g2(xaxis, yaxis, mng_mr, test::point3{0, 0, 0});
    copy(g2_buffer._buffer, g2.data());

    // Check if each bin has 100 points
    EXPECT_EQ(g2.data()[0].size(), 100);
    EXPECT_EQ(g2.data()[1].size(), 100);
    EXPECT_EQ(g2.data()[2].size(), 100);
    EXPECT_EQ(g2.data()[3].size(), 100);
}
