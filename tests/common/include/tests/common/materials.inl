/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

/// Detray include(s)
#include "detray/materials/homogeneous_surface_material.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/mixture.hpp"
#include "detray/materials/predefined_materials.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using point2 = __plugin::point2<scalar>;
using point3 = __plugin::point3<scalar>;

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

    // @note molar density is obtained by the following equation:
    // molar_density = mass_density [GeV/c²] / molar_mass [GeV/c²] * mol / cm³
    EXPECT_FLOAT_EQ(beryllium_tml<scalar>().molar_density(),
                    1.848 / beryllium_tml<scalar>().Ar() * unit_constants::mol /
                        unit_constants::cm3);

    // silicon
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().X0(), 95.7 * unit_constants::mm);
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().L0(), 465.2 * unit_constants::mm);
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().Ar(), 28.03);
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().Z(), 14.);
    EXPECT_FLOAT_EQ(silicon_tml<scalar>().molar_density(),
                    2.32 / silicon_tml<scalar>().Ar() * unit_constants::mol /
                        unit_constants::cm3);
}

TEST(materials, mixture) {

    // Check if material property doesn't change after mixing with other
    // material of 0 ratio
    using mat1 = oxygen_gas<scalar>;
    using mat2 = mixture<scalar, oxygen_gas<scalar>,
                         aluminium<scalar, std::ratio<0, 1>>>;

    EXPECT_FLOAT_EQ(mat1().X0(), mat2().X0());
    EXPECT_FLOAT_EQ(mat1().L0(), mat2().L0());
    EXPECT_FLOAT_EQ(mat1().Ar(), mat2().Ar());
    EXPECT_FLOAT_EQ(mat1().Z(), mat2().Z());
    EXPECT_FLOAT_EQ(mat1().mass_density(), mat2().mass_density());
    EXPECT_FLOAT_EQ(mat1().molar_density(), mat2().molar_density());

    // Air mixture check
    mixture<scalar, carbon_gas<scalar, std::ratio<0, 100>>,
            nitrogen_gas<scalar, std::ratio<76, 100>>,
            oxygen_gas<scalar, std::ratio<23, 100>>,
            argon_gas<scalar, std::ratio<1, 100>>>
        air_mixture;

    EXPECT_TRUE(std::abs(air_mixture.X0() - air<scalar>().X0()) /
                    air<scalar>().X0() <
                0.01);

    material_slab<material<scalar>> slab1(air_mixture, 5.5);
    material_slab<material<scalar>> slab2(air<scalar>(), 2.3);
    material_slab<material<scalar>> slab3(oxygen_gas<scalar>(), 2);

    homogeneous_surface_material surf1(slab1);
    homogeneous_surface_material surf2(slab2);
    homogeneous_surface_material surf3(slab3);

    // Vector check
    std::vector<homogeneous_surface_material<material_slab<material<scalar>>>>
        surface_vec;
    surface_vec.push_back(surf1);
    surface_vec.push_back(surf2);
    surface_vec.push_back(surf3);

    EXPECT_FLOAT_EQ(surface_vec[0].get(point2({0, 0})).thickness_in_X0(),
                    slab1.thickness() / slab1.material().X0());
    EXPECT_FLOAT_EQ(surface_vec[1].get(point2({0, 0})).thickness_in_X0(),
                    slab2.thickness() / slab2.material().X0());
    EXPECT_FLOAT_EQ(surface_vec[2].get(point2({0, 0})).thickness_in_X0(),
                    slab3.thickness() / slab3.material().X0());
}
