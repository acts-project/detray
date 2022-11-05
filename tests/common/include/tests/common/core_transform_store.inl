/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/detail/single_store.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a static transform store
TEST(ALGEBRA_PLUGIN, static_transform_store) {
    using namespace detray;
    using transform3 = __plugin::transform3<detray::scalar>;
    using point3 = __plugin::point3<detray::scalar>;

    using transform_store_t = single_store<transform3>;
    transform_store_t static_store;
    typename transform_store_t::context_type ctx0{};
    typename transform_store_t::context_type ctx1{};

    ASSERT_TRUE(static_store.empty(ctx0));

    ASSERT_EQ(static_store.size(ctx0), 0u);

    point3 t0{0.f, 0.f, 0.f};
    transform3 tf0{t0};
    static_store.push_back(tf0, ctx0);
    ASSERT_EQ(static_store.size(ctx0), 1u);

    point3 t1{1.f, 0.f, 0.f};
    transform3 tf1{t1};
    static_store.push_back(tf1, ctx1);
    ASSERT_EQ(static_store.size(ctx1), 2u);

    point3 t2{2.f, 0.f, 0.f};
    transform3 tf2{t2};
    static_store.push_back(std::move(tf2), ctx0);
    ASSERT_EQ(static_store.size(ctx0), 3u);

    point3 t3{2.f, 0.f, 0.f};
    static_store.emplace_back(ctx0, std::move(t3));
    ASSERT_EQ(static_store.size(ctx0), 4u);

    static_store.emplace_back(ctx0);
    ASSERT_EQ(static_store.size(ctx0), 5u);
}
