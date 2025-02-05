# SPDX-PackageName: "detray, a part of the ACTS project"
# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: MPL-2.0

# Set up the covfie tests.
detray_add_unit_test( covfie
   "constant_field.cpp"
   LINK_LIBRARIES GTest::gtest_main covfie::core
)
