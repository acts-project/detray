/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#pragma once

#define EXPECT_POINT3_NEAR(point1, point2, abs_error) \
    do {                                              \
        EXPECT_NEAR(point1[0], point2[0], abs_error); \
        EXPECT_NEAR(point1[1], point2[1], abs_error); \
        EXPECT_NEAR(point1[2], point2[2], abs_error); \
    } while (false)
