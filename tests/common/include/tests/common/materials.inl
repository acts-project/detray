/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

/// Detray include(s)
#include "detray/materials/material.hpp"
#include "detray/materials/material_composition.hpp"
#include "detray/materials/material_list.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

// This tests the material functionalities

TEST(materials, material) {
    // vacuum
    EXPECT_FLOAT_EQ(vacuum<scalar>().X0(),
                    std::numeric_limits<scalar>::infinity());
    EXPECT_FLOAT_EQ(vacuum<scalar>().L0(),
                    std::numeric_limits<scalar>::infinity());
    EXPECT_FLOAT_EQ(vacuum<scalar>().Ar(), 0);
    EXPECT_FLOAT_EQ(vacuum<scalar>().Z(), 0);
    EXPECT_FLOAT_EQ(vacuum<scalar>().molar_density(), 0);
    EXPECT_FLOAT_EQ(vacuum<scalar>().molar_electron_density(), 0);

    // beryllium
    EXPECT_FLOAT_EQ(beryllium_tml<scalar>().X0(), 352.8 * unit_constants::mm);
    EXPECT_FLOAT_EQ(beryllium_tml<scalar>().L0(), 407.0 * unit_constants::mm);
    EXPECT_FLOAT_EQ(beryllium_tml<scalar>().Ar(), 9.012);
    EXPECT_FLOAT_EQ(beryllium_tml<scalar>().Z(), 4.0);
    EXPECT_FLOAT_EQ(beryllium_tml<scalar>().molar_density(),
                    1.848 / beryllium_tml<scalar>().Ar() * unit_constants::mol /
                        unit_constants::cm3);

    // silicon
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().X0(), 95.7 * unit_constants::mm);
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().L0(), 465.2 * unit_constants::mm);
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().Ar(), 28.03);
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().Z(), 14.);

    // @note molar density is obtained by the following equation:
    // molar_density = 2.32 [GeV/c^2] / molar_mass [GeV/c^2] * mol / cm^3
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().molar_density(),
                    2.32 / silicon_tml<scalar>().Ar() * unit_constants::mol /
                        unit_constants::cm3);
}

TEST(materials, material_composition) {

    material<scalar, std::ratio<1, 6>> m1(0, 0, 0, 0, 0);
    material<scalar, std::ratio<2, 6>> m2(0, 0, 0, 0, 0);
    material<scalar, std::ratio<1, 2>> m3(0, 0, 0, 0, 0);

    material_composition c1(m1, m2, m3);
}