/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <cstdlib>
#include <iostream>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

#include "surfaces_finder_cuda_kernel.hpp"

TEST(surfaces_finder_cuda, surfaces_finder) {

    /** cuda managed memory resource **/
    vecmem::cuda::managed_memory_resource mng_mr;

    /** Make surface finder object **/
    surfaces_finder<n_grids, darray, thrust::tuple, vecmem::vector,
                    vecmem::jagged_vector>
        finder(mng_mr);

    finder[0] = surfaces_regular_circular_grid_t(
        axis::regular<>{2, 0., 2., mng_mr}, axis::circular<>{3, 2., 4., mng_mr},
        mng_mr);

    finder[1] = surfaces_regular_circular_grid_t(
        axis::regular<>{2, 1., 3., mng_mr}, axis::circular<>{2, 2., 4., mng_mr},
        mng_mr);

    /** Fill the grids **/
    for (unsigned int i_g = 0; i_g < finder.size(); i_g++) {
        auto& g2 = finder[i_g];
        auto& xaxis = g2.axis_p0();
        auto& yaxis = g2.axis_p1();

        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
            for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {
                g2.populate(i_x, i_y, i_x + i_y * xaxis.bins());
            }
        }
    }

    /** Output vector to read grid elements **/
    vecmem::vector<dindex> outputs(10, &mng_mr);

    /** Get the data object **/
    auto finder_data = get_data(finder, mng_mr);
    auto outputs_data = vecmem::get_data(outputs);

    /** Run the test kernel **/
    surfaces_finder_test(finder_data, outputs_data);

    /** Check the values **/
    EXPECT_EQ(outputs[0], 0u);
    EXPECT_EQ(outputs[1], 1u);
    EXPECT_EQ(outputs[2], 2u);
    EXPECT_EQ(outputs[3], 3u);
    EXPECT_EQ(outputs[4], 4u);
    EXPECT_EQ(outputs[5], 5u);
    EXPECT_EQ(outputs[6], 0u);
    EXPECT_EQ(outputs[7], 1u);
    EXPECT_EQ(outputs[8], 2u);
    EXPECT_EQ(outputs[9], 3u);
}
