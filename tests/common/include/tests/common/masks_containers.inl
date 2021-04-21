/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/rectangle2.hpp"
#include "masks/trapezoid2.hpp"
#include "masks/cylinder3.hpp"
#include "masks/single3.hpp"
#include "utils/containers.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace __plugin;

// This tests the construction of a surface
TEST(mask, tuple)
{
    using rectangle_masks = dvector<rectangle2<>>;
    using cylinder_masks = dvector<cylinder3<>>;

    rectangle_masks rectangles;
    cylinder_masks cylinders;
    dtuple< rectangle_masks, cylinder_masks > all_masks;

    std::get<rectangle_masks>(all_masks).push_back( {3.,4.});
    std::get<rectangle_masks>(all_masks).push_back( {2.,8.});

    auto intersector = std::get<rectangle_masks>(all_masks)[0].intersector();
    
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}