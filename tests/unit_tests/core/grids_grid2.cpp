/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "utils/indexing.hpp"

#include <gtest/gtest.h>

#include <climits>

using namespace detray;

TEST(grids, grid2_replace_populator)
{
    replace_populator<guaranteed_index, std::numeric_limits<guaranteed_index>::max()> replacer;
    serializer2 serializer;

    axis::closed<> xaxis{10, -5., 5.};
    axis::closed<> yaxis{10, -5., 5.};
    using grid2r = grid2<decltype(replacer), decltype(xaxis), decltype(yaxis), decltype(serializer)>;

    grid2r g2(std::move(xaxis), std::move(yaxis));

    // Test the initialization
    test::point2 p = {-4.5, -4.5};
    for (unsigned int ib0 = 0; ib0 < 10; ++ib0)
    {
        for (unsigned int ib1 = 0; ib1 < 10; ++ib1)
        {
            p = {-4.5 + ib0, -4.5 + ib1};
            EXPECT_EQ(g2.bin(p), std::numeric_limits<guaranteed_index>::max());
        }
    }

    p = {-4.5, -4.5};
    // Fill and read
    g2.populate(p, 3u);
    EXPECT_EQ(g2.bin(p), 3u);

    // Fill and read two times, fill first 0-99, then 100-199
    for (unsigned int il = 0; il < 2; ++il)
    {
        unsigned int counter = il * 100;
        for (unsigned int ib0 = 0; ib0 < 10; ++ib0)
        {
            for (unsigned int ib1 = 0; ib1 < 10; ++ib1)
            {
                p = {-4.5 + ib0, -4.5 + ib1};
                g2.populate(p, counter);
                EXPECT_EQ(g2.bin(p), counter++);
            }
        }
    }

    // A zone test w/o neighbour hood 
    p = {-4.5, -4.5};
    auto test = g2.zone(p, {0, 0}, true);
    dvector<guaranteed_index> expect = { 100u };
    EXPECT_EQ(test, expect);


    // A zone test with neighbour hood
    p = {0.5, 0.5};
    test = g2.zone(p, {1, 2}, true);
    expect = {143u, 144u, 145u, 146u, 147u, 153u, 154u, 155u, 156u, 157u, 163u, 164u, 165u, 166u, 167u};
    EXPECT_EQ(test, expect);

    axis::circular<> circular{4, -2., 2.};
    axis::closed<> closed{5, 0., 5.};
    using grid2cc = grid2<decltype(replacer), decltype(circular), decltype(closed), decltype(serializer)>;

    grid2cc g2cc(std::move(circular), std::move(closed));
    unsigned int counter = 0;
    for (unsigned icl = 0; icl < 5; ++icl)
    {
        for (unsigned ici = 0; ici < 4; ++ici)
        {
            p = {-1.5 + ici, 0.5 + icl};
            g2cc.populate(p, counter++);
        }
    }

    // A zone test for circular testing
    p = {1.5, 2.5};
    test = g2cc.zone(p, {1, 1}, true);
    expect = {4u, 6u, 7u, 8u, 10u, 11u, 12u, 14u, 15u};
    EXPECT_EQ(test, expect);
}

TEST(grids, grid2_complete_populator)
{

    complete_populator<optional_index, -1, 3, false> completer;
    serializer2 serializer;

    axis::closed<> xaxis{2, -1., 1.};
    axis::closed<> yaxis{2, -1., 1.};
    using grid2r = grid2<decltype(completer), decltype(xaxis), decltype(yaxis), decltype(serializer)>;

    grid2r g2(std::move(xaxis), std::move(yaxis));

    // Test the initialization
    test::point2 p = {-0.5, -0.5};
    decltype(completer)::store_value invalid = {-1, -1, -1};
    for (unsigned int ib0 = 0; ib0 < 2; ++ib0)
    {
        for (unsigned int ib1 = 0; ib1 < 2; ++ib1)
        {
            p = {-0.5 + ib0, -0.5 + ib1};
            EXPECT_EQ(g2.bin(p), invalid);
        }
    }

    // Fill and read
    p = {-0.5, -0.5};
    g2.populate(p, 4u);

    decltype(completer)::store_value expected = {4u, -1, -1};
    auto test = g2.bin(p);
    EXPECT_EQ(test, expected);

    auto zone_test = g2.zone(p, {0, 0});
    dvector<optional_index> zone_expected = {4u};
    EXPECT_EQ(zone_test, zone_expected);

    g2.populate(p, 2u);
    expected = {4u, 2u, -1};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    g2.populate(p, 7u);
    expected = {4u, 2u, 7u};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    // Bin is completed, new entry is ignored
    g2.populate(p, 16u);
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    // Zone test of a complete bin
    zone_test = g2.zone(p, {0, 0});
    zone_expected = {4u, 2u, 7u};
    EXPECT_EQ(zone_test, zone_expected);

    // Fill some other bins
    p = {0.5, -0.5};
    g2.populate(p, 16u);

    p = {0.5, 0.5};
    g2.populate(p, 17u);
    g2.populate(p, 18u);

    zone_test = g2.zone(p, {1, 1});
    zone_expected = {4u, 2u, 7u, 16u, 17u, 18u};
    EXPECT_EQ(zone_test, zone_expected);
}

TEST(grids, grid2_attach_populator)
{

    attach_populator<guaranteed_index> attacher;
    serializer2 serializer;

    axis::closed<> xaxis{2, -1., 1.};
    axis::closed<> yaxis{2, -1., 1.};
    using grid2r = grid2<decltype(attacher), decltype(xaxis), decltype(yaxis), decltype(serializer)>;

    grid2r g2(std::move(xaxis), std::move(yaxis));

    // Test the initialization
    test::point2 p = {-0.5, -0.5};
    decltype(attacher)::store_value invalid = {};
    for (unsigned int ib0 = 0; ib0 < 2; ++ib0)
    {
        for (unsigned int ib1 = 0; ib1 < 2; ++ib1)
        {
            p = {-0.5 + ib0, -0.5 + ib1};
            EXPECT_EQ(g2.bin(p), invalid);
        }
    }

    p = {-0.5, -0.5};
    g2.populate(p, 4u);

    decltype(attacher)::store_value expected = {4u};
    auto test = g2.bin(p);
    EXPECT_EQ(test, expected);

    auto zone_test = g2.zone(p, {0, 0});
    dvector<guaranteed_index> zone_expected = {4u};
    EXPECT_EQ(zone_test, zone_expected);

    p = {-0.5, 0.5};
    g2.populate(p, 9u);

    p = {0.5, -0.5};
    g2.populate(p, 1u);

    p = {0.5, 0.5};
    g2.populate(p, 7u);

    expected = {7u};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    zone_test = g2.zone(p, {1, 1}, true);
    zone_expected = {1u, 4u, 7u, 9u};
    EXPECT_EQ(zone_test, zone_expected);
}


TEST(grids, grid2_shift)
{
    replace_populator<guaranteed_index, 0> replacer;
    serializer2 serializer;

    axis::closed<> xaxis{10, -5., 5.};
    axis::closed<> yaxis{10, -5., 5.};

    using grid2r = grid2<decltype(replacer), decltype(xaxis), decltype(yaxis), decltype(serializer)>;

    grid2r g2(std::move(xaxis), std::move(yaxis));

    // Test the initialization
    test::point2 p = {-4.5, -4.5};
    EXPECT_EQ(g2.bin(p), 0u);

    g2.shift(8u);
    EXPECT_EQ(g2.bin(p), 8u);

}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
