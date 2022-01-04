/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>

// detray test
#include "tests/common/test_defs.hpp"
#include "tests/common/tools/hash_tree.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the correctness of the hash tree tool
TEST(ALGEBRA_PLUGIN, hash_tree) {
    using namespace detray;

    dvector<dindex> test_matrix_1 = {1, 1, 1, 1, 1, 1};

    auto ht_1 = hash_tree(test_matrix_1);
    decltype(ht_1)::hash_function hasher{};
    using hash_type = decltype(ht_1)::hash_type;

    dvector<hash_type> digests{};
    for (auto input : test_matrix_1) {
        digests.push_back(hasher(input));
    }
    digests.push_back(0);


    std::cout << ht_1.to_string() << std::endl;
    
    dvector<dindex> adj_mat_truth = {2, 3, 0, 4, 2, 0, 0, 0, 0};

    // Check this with graph
    ASSERT_EQ(ht_1.root(), 6);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
