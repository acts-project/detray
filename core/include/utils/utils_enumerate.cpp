/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "utils/enumerate.hpp"

#include <gtest/gtest.h>

using namespace detray;

// This tests the convenience enumeration function
TEST(utils, enumerate_function) {

  struct uint_holder {
    unsigned int ui = 0;
  };

  dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}};

  using container_type_iter =
      decltype(std::begin(std::declval<dvector<uint_holder>>()));

  for (auto [i, v] : enumerate(seq)) {
    ASSERT_EQ(i, v.ui);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}