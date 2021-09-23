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

#include "grids_grid2_cuda_kernel.cuh"

using namespace detray;

TEST(grids_cuda, grid2_replace_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // axis
    axis::regular<> xaxis{4, -1., 3.};
    axis::regular<> yaxis{6, 0., 6.};

    auto x_interval = (xaxis.max - xaxis.min) / xaxis.n_bins;
    auto y_interval = (yaxis.max - yaxis.min) / yaxis.n_bins;

    // declare grid
    grid2r_replace g2(std::move(xaxis), std::move(yaxis), mng_mr,
                      test::point3{0, 0, 0});

    // get grid_data
    grid2_data g2_data(g2);

    // pre-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {

            auto& data = g2.bin(i_x, i_y);

            EXPECT_EQ(data, g2.populator().kInvalid);
        }
    }

    // run test function in cuda
    grid_test1(g2_data);

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
    grid2r_complete g2(std::move(xaxis), std::move(yaxis), mng_mr,
                       test::point3{0, 0, 0});

    // get grid_data
    grid2_data g2_data(g2);

    // pre-check
    for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {

            auto& data = g2.bin(i_x, i_y);

            for (auto pt : data) {
                EXPECT_EQ(pt, g2.populator().kInvalid);
            }
        }
    }

    // run test function in cuda
    grid_test2(g2_data);

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

TEST(grids_cuda, grid2_attach_populator) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // axis
    axis::regular<> xaxis{4, -1., 3.};
    axis::regular<> yaxis{6, 0., 6.};

    auto x_interval = (xaxis.max - xaxis.min) / xaxis.n_bins;
    auto y_interval = (yaxis.max - yaxis.min) / yaxis.n_bins;

    // declare grid
    grid2r_attach g2(std::move(xaxis), std::move(yaxis), mng_mr,
                     test::point3{0, 0, 0});

    // fill the grid
    for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
        for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
            for (int i_p = 0; i_p < n_points; i_p++) {

                auto bin_id = i_x + i_y * xaxis.bins();
                auto gid = i_p + bin_id * n_points;

                test::point3 tp({xaxis.min + gid * x_interval,
                                 yaxis.min + gid * y_interval, 0.5});
                g2.populate(i_x, i_y, std::move(tp));
            }
        }
    }

    grid2r_attach::populator_t::serialized_storage original_storage = g2.data();

    // get grid_data
    grid2_data g2_data(g2, &mng_mr);

    // run test function in cuda
    grid_test2(g2_data);

    // check if the vector remains the same
    for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
        for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
            auto bin_id = i_x + i_y * xaxis.bins();
            auto& original_data = original_storage[bin_id];

            auto& data = g2.bin(i_x, i_y);
            for (int i_p = 0; i_p < data.size(); i_p++) {
                auto& pt = data[i_p];
                auto& original_pt = original_data[i_p];

                EXPECT_EQ(pt, original_pt);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
