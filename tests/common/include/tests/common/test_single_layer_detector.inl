/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "tools/navigator.hpp"
#include "tests/common/single_layer_detector.hpp"

#include <fstream>
#include <cmath>
#include <climits>

#include <gtest/gtest.h>

unsigned int theta_steps = 10;
unsigned int phi_steps = 1000;
bool stream_file = true;

auto d = createDetector();

// This test navigates through a cylindrical detector
TEST(__plugin, construction)
{
    // Let's inspect the detector
    // 2 volumes
    ASSERT_EQ(d.volumes().size(), 2u);

    // 6 portals: 2 bp ec, 1 shared, 2 barrrel ec, 1 barrel cover
    ASSERT_EQ(d.portal_surfaces().size(), 6u);


}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
