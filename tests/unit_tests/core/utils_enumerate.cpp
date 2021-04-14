/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "utils/indexing.hpp"
#include "utils/containers.hpp"
#include "utils/enumerate.hpp"
#include "tests/common/test_defs.hpp"

#include <gtest/gtest.h>


using namespace detray;

// This tests the convenience range_enumeration function: single
TEST(utils, sequence_single)
{

    dindex check = 0;
    dindex single = 7;
    for (auto i : sequence(single)){
        check += i;
    }
    ASSERT_EQ(check, single);

}

// This tests the convenience range_enumeration function: range
TEST(utils, sequence_range)
{

    darray<dindex, 2> range = {2,7};
    std::vector<dindex> reference = {2, 3, 4, 5, 6, 7};
    std::vector<dindex> check = {};
    for (auto i : sequence(range)){
        check.push_back(i);
    }
    ASSERT_EQ(check, reference);

}

// This tests the convenience enumeration function
TEST(utils, enumerate)
{

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = { {0}, {1}, {2}, {3}, {4}, {5}};

    using container_type_iter = decltype(std::begin(std::declval<dvector<uint_holder>>()));

    for (auto [i, v] : enumerate(seq)){
        ASSERT_EQ(i, v.ui);
    }

}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}