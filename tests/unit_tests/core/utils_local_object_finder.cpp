/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <functional>
#include <vecmem/memory/host_memory_resource.hpp>

#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/populator.hpp"
#include "grids/serializer2.hpp"
#include "tests/common/test_defs.hpp"
#include "tools/local_object_finder.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

using namespace detray;

// This tests the convenience enumeration function
TEST(utils, local_object_finder) {
    vecmem::host_memory_resource host_mr;

    replace_populator<> replacer;
    serializer2 serializer;

    test::point2 p2 = {-4.5, -4.5};

    axis::regular<> xaxis{10, -5., 5.};
    axis::regular<> yaxis{10, -5., 5.};
    using grid2r = grid2<decltype(replacer), decltype(xaxis), decltype(yaxis),
                         decltype(serializer)>;

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr);

    g2.populate(p2, 8u);

    dvector<dindex> expected = {8u};
    EXPECT_EQ(g2.zone(p2), expected);

    local_zone_finder<grid2r> local_zone(std::move(g2));

    local_single_finder<dindex> local_single(8u);

    using local_finder = std::function<dvector<dindex>(const test::point2 &)>;

    std::vector<local_finder> local_finders = {local_zone, local_single};

    EXPECT_EQ(local_finders[0](p2), expected);
    EXPECT_EQ(local_finders[1](p2), expected);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
