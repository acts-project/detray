/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "grids/populator.hpp"
#include "utils/indexing.hpp"

#include <gtest/gtest.h>

#include <climits>

using namespace detray;

TEST(grids, replace_populator)
{
    replace_populator<guaranteed_index> replacer;
    guaranteed_index stored = 3;
    replacer(stored, 2);
    EXPECT_EQ(stored, 2u);

    replacer(stored, 42);
    EXPECT_EQ(stored, 42u);
}

TEST(grids, complete_populator)
{

    using cpopulator4 = complete_populator<guaranteed_index, std::numeric_limits<guaranteed_index>::max(), 4>;
    cpopulator4 completer;

    cpopulator4::store_value stored = {cpopulator4::invalid_value, cpopulator4::invalid_value, cpopulator4::invalid_value, cpopulator4::invalid_value};


    cpopulator4::store_value test = stored;
    test[0] = 9u;
    completer(stored, 9);
    EXPECT_EQ(stored, test);

    test[1] = 3u;
    completer(stored, 3);
    EXPECT_EQ(stored, test);

    using sort_cpopulator4 = complete_populator<guaranteed_index, std::numeric_limits<guaranteed_index>::max(), 4, true>;
    sort_cpopulator4 sort_completer;

    test = { 0, 3, 9, 1000};
    sort_completer(stored, 1000);
    sort_completer(stored, 0);
    EXPECT_EQ(stored, test);

}

TEST(grids, attach_populator)
{
    // Attch populator without sorting
    attach_populator<guaranteed_index> attacher;
    attach_populator<guaranteed_index>::store_value stored = {3};
    attacher(stored, 2);
    attach_populator<guaranteed_index>::store_value test = {3u, 2u};
    EXPECT_EQ(stored, test);

    attacher(stored, 42);
    test = {3u, 2u, 42u};
    EXPECT_EQ(stored, test);

    // Attach populator with sorting
    attach_populator<guaranteed_index, true> sort_attacher;
    sort_attacher(stored, 11);
    test = {2u, 3u, 11u, 42u};
    EXPECT_EQ(stored, test);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
