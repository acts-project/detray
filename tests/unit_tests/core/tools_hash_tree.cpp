/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <cmath>

// detray test
#include "tests/common/tools/hash_tree.hpp"

namespace {
template <typename hasher_t, typename digest_collection_t>
void test_hash(std::size_t first_child, std::size_t n_prev_level,
               digest_collection_t &digests) {
    // base case
    if (n_prev_level <= 1u) {
        return;
    }

    auto last_child = first_child + n_prev_level;

    // Run over previous tree level to build the next level
    for (std::size_t i = first_child; i < last_child; i += 2u) {
        auto digest = hasher_t{}(digests[i], digests[i + 1u]);
        digests.push_back(digest);
    }
    std::size_t n_level =
        static_cast<std::size_t>(0.5f * static_cast<float>(n_prev_level));
    // Need dummy leaf node for next level?
    if (n_level % 2u != 0u and n_level > 1u) {
        digests.push_back(0);
        n_level++;
    }
    // begin next time where we ended this time
    test_hash<hasher_t>(last_child, n_level, digests);
}

}  // namespace

// This tests the correctness of the hash tree tool
TEST(ALGEBRA_PLUGIN, hash_tree) {
    using namespace detray;

    dvector<dindex> test_matrix = {1u, 1u, 1u, 1u, 1u, 1u};

    auto ht = hash_tree(test_matrix);
    using hash_tree_t = decltype(ht);
    hash_tree_t::hash_function hasher{};

    dvector<hash_tree_t::hash_type> digests{};
    // Build up first level
    for (auto input : test_matrix) {
        digests.push_back(hasher(input));
    }
    // Build the other levels
    test_hash<hash_tree_t::hash_function>(0, digests.size(), digests);

    // Check this with graph
    EXPECT_EQ(ht.root(), 6);

    const auto tree_data = ht.tree();
    EXPECT_EQ(digests.size(), tree_data.size());
    for (std::size_t n_idx = 0u; n_idx < tree_data.size(); ++n_idx) {
        EXPECT_EQ(digests[n_idx], tree_data[n_idx].key());
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
