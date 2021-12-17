/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/surfaces_finder.hpp"

using namespace detray;

// This test checks the surfaces finder
TEST(ALGEBRA_PLUGIN, surfaces_finder) {

    /** Host memory resource **/
    vecmem::host_memory_resource host_mr;

    /** Make surfaces finder **/
    constexpr int n_grids = 2;

    /** Make surface finder object **/
    surfaces_finder<n_grids> finder(host_mr);
    using surfaces_regular_circular_grid_t =
        decltype(finder)::surfaces_regular_circular_grid;

    finder[0] = surfaces_regular_circular_grid_t(
        axis::regular<>{2, 0., 2., host_mr},
        axis::circular<>{3, 2., 4., host_mr}, host_mr);

    finder[1] = surfaces_regular_circular_grid_t(
        axis::regular<>{2, 1., 5., host_mr},
        axis::circular<>{2, 2., 4., host_mr}, host_mr);

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

    /** Check the values **/
    EXPECT_EQ(finder[0].bin(0, 0)[0], 0u);
    EXPECT_EQ(finder[0].bin(1, 0)[0], 1u);
    EXPECT_EQ(finder[0].bin(0, 1)[0], 2u);
    EXPECT_EQ(finder[0].bin(1, 1)[0], 3u);
    EXPECT_EQ(finder[0].bin(0, 2)[0], 4u);
    EXPECT_EQ(finder[0].bin(1, 2)[0], 5u);

    EXPECT_EQ(finder[1].bin(0, 0)[0], 0u);
    EXPECT_EQ(finder[1].bin(1, 0)[0], 1u);
    EXPECT_EQ(finder[1].bin(0, 1)[0], 2u);
    EXPECT_EQ(finder[1].bin(1, 1)[0], 3u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}