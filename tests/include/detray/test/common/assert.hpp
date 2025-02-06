// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#include <gtest/gtest.h>

#pragma once

#define EXPECT_POINT3_NEAR(point1, point2, abs_error) \
    do {                                              \
        EXPECT_NEAR(point1[0], point2[0], abs_error); \
        EXPECT_NEAR(point1[1], point2[1], abs_error); \
        EXPECT_NEAR(point1[2], point2[2], abs_error); \
    } while (false)
