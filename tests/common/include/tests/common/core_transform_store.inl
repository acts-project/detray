/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/transform_store.hpp"

#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a static transform store
TEST(__plugin, static_transform_store)
{
    using namespace detray;
    using namespace __plugin;

    static_transform_store static_store;   
    static_transform_store::context ctx0;
    static_transform_store::context ctx1;

    ASSERT_TRUE(static_store.empty(ctx0));

    ASSERT_EQ(static_store.size(ctx0), 0u);

    transform3::point3 t0{0.,0.,0.};
    transform3 tf0{t0};
    static_store.push_back(ctx0, tf0);
    ASSERT_EQ(static_store.size(ctx0), 1u);

    transform3::point3 t1{1.,0.,0.};
    transform3 tf1{t1};
    static_store.push_back(ctx1, tf1);
    ASSERT_EQ(static_store.size(ctx1), 2u);

    transform3::point3 t2{2.,0.,0.};
    transform3 tf2{t2};
    static_store.push_back(ctx0, std::move(tf2));
    ASSERT_EQ(static_store.size(ctx0), 3u);

}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

