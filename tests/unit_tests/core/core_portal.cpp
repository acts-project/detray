/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "core/portal.hpp"
#include "core/intersection.hpp"
#include "utils/containers.hpp"

#include <exception>
#include <gtest/gtest.h>

using namespace detray;
using namespace plugin;

namespace
{
    struct dummy_mask
    {
        point3 _v = {0., 0., 0.};

        const intersection_status operator()(const point3 &p) const
        {
            if (_v == p)
            {
                return intersection_status::e_inside;
            }
            return intersection_status::e_outside;
        }
    };
} // namespace

// This tests a portal with one link each side
TEST(core, portal_single_mask)
{
    dummy_mask along_mask = {1., 2., 3.};
    dummy_mask opposite_mask = {4., 5., 6.};

    using single_mask = darray<dummy_mask, 1>;
    using single_offset = darray<int, 1>;

    single_mask along_container = {along_mask};
    single_offset along_offset = {7};
    single_mask opposite_container = {opposite_mask};
    single_offset opposite_offset = {42};

    portal<unsigned int, single_mask, single_offset> single_portal(
        0, std::move(along_container), std::move(along_offset), std::move(opposite_container), std::move(opposite_offset));

    point3 along_hit = {1., 2., 3.};
    point3 opposite_hit = {4., 5., 6.};

    ASSERT_EQ(single_portal.surface(), 0);

    ASSERT_EQ(single_portal.along(along_hit), 7);
    EXPECT_THROW(single_portal.along(opposite_hit), std::invalid_argument);

    EXPECT_THROW(single_portal.opposite(along_hit), std::invalid_argument);
    ASSERT_EQ(single_portal.opposite(opposite_hit), 42);
}

// This tests a portal with different numbers of links each side
TEST(core, portal_multi_mask)
{
    dummy_mask along_mask = {1., 2., 3.};
    dummy_mask opposite_mask = {4., 5., 6.};

    using mask_container = dvector<dummy_mask>;
    using offset_container = dvector<int>;

    dummy_mask a0, a2, a3, o0, o1;

    mask_container along_container = {a0, along_mask, a2, a3};
    offset_container along_offset = {3, 7, 2, 1};
    mask_container opposite_container = {o0, o1, opposite_mask};
    offset_container opposite_offset = {-1, -8, 42};

    portal<unsigned int, mask_container, offset_container> multi_portal(
        1, std::move(along_container), std::move(along_offset), std::move(opposite_container), std::move(opposite_offset));

    point3 along_hit = {1., 2., 3.};
    point3 opposite_hit = {4., 5., 6.};

    ASSERT_EQ(multi_portal.surface(), 1);

    ASSERT_EQ(multi_portal.along(along_hit), 8); // includes offset
    EXPECT_THROW(multi_portal.along(opposite_hit), std::invalid_argument);

    EXPECT_THROW(multi_portal.opposite(along_hit), std::invalid_argument);
    ASSERT_EQ(multi_portal.opposite(opposite_hit), 44); // includes offset
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
