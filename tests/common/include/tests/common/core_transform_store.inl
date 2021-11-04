/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/transform_store.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a static transform store
TEST(ALGEBRA_PLUGIN, static_transform_store) {
    using namespace detray;
    using namespace __plugin;

    static_transform_store<> static_store;
    static_transform_store<>::context ctx0;
    static_transform_store<>::context ctx1;

    ASSERT_TRUE(static_store.empty(ctx0));

    ASSERT_EQ(static_store.size(ctx0), 0u);

    point3 t0{0., 0., 0.};
    transform3 tf0{t0};
    static_store.push_back(ctx0, tf0);
    ASSERT_EQ(static_store.size(ctx0), 1u);

    point3 t1{1., 0., 0.};
    transform3 tf1{t1};
    static_store.push_back(ctx1, tf1);
    ASSERT_EQ(static_store.size(ctx1), 2u);

    point3 t2{2., 0., 0.};
    transform3 tf2{t2};
    static_store.push_back(ctx0, std::move(tf2));
    ASSERT_EQ(static_store.size(ctx0), 3u);

    point3 t3{2., 0., 0.};
    static_store.emplace_back(ctx0, std::move(t3));
    ASSERT_EQ(static_store.size(ctx0), 4u);

    static_store.emplace_back(ctx0);
    ASSERT_EQ(static_store.size(ctx0), 5u);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
