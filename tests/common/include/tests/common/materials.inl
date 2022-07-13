/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

/// Detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/material_rod.hpp"
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
        air_mix;

    EXPECT_TRUE(std::abs(air_mix.X0() - air<scalar>().X0()) /
                    air<scalar>().X0() <
                0.01);
    EXPECT_TRUE(std::abs(air_mix.L0() - air<scalar>().L0()) /
                    air<scalar>().L0() <
                0.01);
    EXPECT_TRUE(std::abs(air_mix.Ar() - air<scalar>().Ar()) /
                    air<scalar>().Ar() <
                0.01);
    EXPECT_TRUE(std::abs(air_mix.Z() - air<scalar>().Z()) / air<scalar>().Z() <
                0.01);
    EXPECT_TRUE(
        std::abs(air_mix.mass_density() - air<scalar>().mass_density()) /
            air<scalar>().mass_density() <
        0.01);

    // Vector check
    material_slab<scalar> slab1(air_mix, 5.5);
    material_slab<scalar> slab2(air<scalar>(), 2.3);
    material_slab<scalar> slab3(oxygen_gas<scalar>(), 2);

    std::vector<material_slab<scalar>> slab_vec;

    slab_vec.push_back(slab1);
    slab_vec.push_back(slab2);
    slab_vec.push_back(slab3);

    EXPECT_FLOAT_EQ(slab_vec[0].thickness_in_X0(),
                    slab1.thickness() / slab1.get_material().X0());
    EXPECT_FLOAT_EQ(slab_vec[1].thickness_in_X0(),
                    slab2.thickness() / slab2.get_material().X0());
    EXPECT_FLOAT_EQ(slab_vec[2].thickness_in_X0(),
                    slab3.thickness() / slab3.get_material().X0());
}

// This tests the material slab functionalities
TEST(materials, material_slab) {

    material_slab<scalar> slab(oxygen_gas<scalar>(),
                               scalar(2) * scalar(unit_constants::mm));

    line_plane_intersection is;
    is.cos_incidence_angle = scalar(0.3);

    EXPECT_FLOAT_EQ(slab.path_segment(is),
                    scalar(2) * scalar(unit_constants::mm) / scalar(0.3));
    EXPECT_FLOAT_EQ(slab.path_segment_in_X0(is),
                    slab.path_segment(is) / slab.get_material().X0());
    EXPECT_FLOAT_EQ(slab.path_segment_in_L0(is),
                    slab.path_segment(is) / slab.get_material().L0());
}

// This tests the material rod functionalities
TEST(materials, material_rod) {

    material_rod<scalar> rod(oxygen_gas<scalar>(),
                             scalar(2) * scalar(unit_constants::mm));

    line_plane_intersection is;
    is.p2[0] = 1. * unit_constants::mm;

    EXPECT_FLOAT_EQ(rod.path_segment(is), scalar(2.) * scalar(std::sqrt(3)));
    EXPECT_FLOAT_EQ(rod.path_segment_in_X0(is),
                    rod.path_segment(is) / rod.get_material().X0());
    EXPECT_FLOAT_EQ(rod.path_segment_in_L0(is),
                    rod.path_segment(is) / rod.get_material().L0());
}