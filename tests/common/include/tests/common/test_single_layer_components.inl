/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/populator.hpp"
#include "grids/serializer2.hpp"
#include "tools/concentric_cylinder_intersector.hpp"
#include "tools/planar_intersector.hpp"
#include "tools/navigator.hpp"
#include "tests/common/single_layer_detector.hpp"
#include "tests/common/test_surfaces.hpp"

#include <fstream>
#include <cmath>
#include <climits>
#include <random>

#include <gtest/gtest.h>

bool write_files = false;
bool test_against_all = true;
unsigned int tests = 100000;

scalar r = 32.;
scalar vr_inner = 29.;
scalar vr_outer = 36.;
scalar vz_half = 400.;
scalar r_stag = 0.25;
unsigned int n_phi = 12;
scalar t_phi = 0.12;
scalar ov_rphi = 0.25;
scalar z_length = 600.;
scalar z_overlap = 2.;
unsigned int n_z = 7;

// Retrieving the barrel components for testing
auto barrel_components = detray::create_barrel_components(r, r_stag, n_phi,
                                                          t_phi, ov_rphi, z_length, z_overlap, n_z,
                                                          vr_inner, vr_outer, vz_half);

auto barrel_module = std::get<0>(barrel_components);
auto barrel_transforms = std::get<1>(barrel_components);
auto barrel_finders = std::get<2>(barrel_components);

// Three-dimensional definitions
using transform3 = __plugin::transform3;
using vector3 = __plugin::transform3::vector3;
using point3 = __plugin::transform3::point3;

// Intersectors
concentric_cylinder_intersector cci;
planar_intersector pi;

// Test the basic barrel components function
TEST(__plugin, barrel_component_construction)
{
    // Rectangular module parameters
    ASSERT_EQ(barrel_module.size(), 2u);
    // Number of Barrel transforms
    ASSERT_EQ(barrel_transforms.size(), n_phi * n_z);
    // Barrel finders
    ASSERT_EQ(barrel_finders.size(), 4u);

    // Check if the finders make sense
    auto inner_finder = barrel_finders[0];
    auto outer_finder = barrel_finders[1];
    auto nec_finder = barrel_finders[2];
    auto pec_finder = barrel_finders[3];

    using cylinder_grid = grid2<replace_populator<>, axis::circular<>, axis::closed<>, serializer2>;
    using disc_grid = grid2<replace_populator<>, axis::closed<>, axis::circular<>, serializer2>;

    using cylinder_finder = local_zone_finder<cylinder_grid>;
    using disc_finder = local_zone_finder<disc_grid>;

    auto inner_local_finder = inner_finder.target<cylinder_finder>();
    auto outer_local_finder = outer_finder.target<cylinder_finder>();
    auto nec_local_finder = nec_finder.target<disc_finder>();
    auto pec_local_finder = pec_finder.target<disc_finder>();

    EXPECT_TRUE(inner_local_finder != nullptr);
    EXPECT_TRUE(outer_local_finder != nullptr);
    EXPECT_TRUE(nec_local_finder != nullptr);
    EXPECT_TRUE(pec_local_finder != nullptr);

    const auto &inner_grid = inner_local_finder->grid();
    auto inner_grid_axis_p0 = inner_grid.axis_p0();
    auto inner_grid_axis_p1 = inner_grid.axis_p1();
    ASSERT_EQ(inner_grid_axis_p0.bins, n_phi);
    ASSERT_EQ(inner_grid_axis_p1.bins, n_z);

    const auto &outer_grid = outer_local_finder->grid();
    auto outer_grid_axis_p0 = outer_grid.axis_p0();
    auto outer_grid_axis_p1 = outer_grid.axis_p1();
    ASSERT_EQ(outer_grid_axis_p0.bins, n_phi);
    ASSERT_EQ(outer_grid_axis_p1.bins, n_z);

    const auto &ndisc_grid = nec_local_finder->grid();
    auto ndisc_grid_axis_p0 = ndisc_grid.axis_p0();
    auto ndisc_grid_axis_p1 = ndisc_grid.axis_p1();
    ASSERT_EQ(ndisc_grid_axis_p0.bins, 1u);
    ASSERT_EQ(ndisc_grid_axis_p1.bins, n_phi);

    const auto &pdisc_grid = pec_local_finder->grid();
    auto pdisc_grid_axis_p0 = pdisc_grid.axis_p0();
    auto pdisc_grid_axis_p1 = pdisc_grid.axis_p1();
    ASSERT_EQ(pdisc_grid_axis_p0.bins, 1u);
    ASSERT_EQ(pdisc_grid_axis_p1.bins, n_phi);

    using cyl_point2 = __plugin::cylindrical2::point2;
    using pol_point2 = __plugin::polar2::point2;

    scalar z_step = 2 * barrel_module[1] - z_overlap;
    scalar phi_step = 2 * M_PI / n_phi;
    unsigned int imod = 0;

    // Test the entries of the inner/outer Barrel grids
    for (unsigned int iz = 0; iz < n_z; ++iz)
    {
        for (unsigned int iphi = 0; iphi < n_phi; ++iphi)
        {

            scalar phi = -M_PI + iphi * phi_step;
            scalar z = -0.5 * z_length + barrel_module[1] + iz * z_step;
            cyl_point2 inner_p2 = {vr_inner * phi, z};
            cyl_point2 outer_p2 = {vr_outer * phi, z};

            auto inner_zone = inner_finder(inner_p2, {1, 1});
            auto outer_zone = outer_finder(outer_p2, {1, 1});

            if (iz == 0)
            {
                ASSERT_EQ(inner_zone.size(), 6u);
                ASSERT_EQ(outer_zone.size(), 6u);
                // nec ring
                pol_point2 nec_p2 = {0.5 * (vr_inner + vr_outer), phi};
                auto nec_zone = nec_finder(nec_p2, {1, 1});
                ASSERT_EQ(nec_zone.size(), 3u);
            }
            else if ((iz + 1) == n_z)
            {
                ASSERT_EQ(inner_zone.size(), 6u);
                ASSERT_EQ(outer_zone.size(), 6u);
            }
            else
            {
                ASSERT_EQ(inner_zone.size(), 9u);
                ASSERT_EQ(outer_zone.size(), 9u);
                // pec ring
                pol_point2 pec_p2 = {0.5 * (vr_inner + vr_outer), phi};
                auto pec_zone = pec_finder(pec_p2, {1, 1});
                ASSERT_EQ(pec_zone.size(), 3u);
            }
            ASSERT_EQ(inner_zone, outer_zone);
        }
    }
}

TEST(__plugin, barrel_object_finding)
{

    std::ofstream search_points;
    std::ofstream heuristic_points;
    std::ofstream confirmed_points;

    if (write_files)
    {
        search_points.open("barrel_layer_search_points.csv");
        search_points << "x,y,z"
                      << "\n";

        confirmed_points.open("barrel_layer_confirmed_points.csv");
        confirmed_points << "x,y,z"
                         << "\n";

        heuristic_points.open("barrel_layer_heuristic_points.csv");
        heuristic_points << "x,y,z"
                         << "\n";
    }

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<scalar> dist_phi(-M_PI, M_PI);
    std::uniform_real_distribution<scalar> dist_theta(0.1, M_PI - 0.1);

    transform3 identity(vector3{0., 0., 0.});
    cylinder3<false> cylinder = {vr_inner, -vz_half, vz_half};

    // Cylindrical and Cartesian local frames
    __plugin::cylindrical2 cylindrical2;
    __plugin::cartesian2 cart2;

    vector3 ori = {0., 0., 0.};

    // Check if the finders make sense
    auto inner_finder = barrel_finders[0];

    rectangle2<> rect = {barrel_module[0], barrel_module[1]};

    // Cross check
    unsigned int hit_modules_found = 0;
    unsigned int hit_modules_heuristic = 0;

    for (unsigned int itest = 0; itest < tests; ++itest)
    {

        scalar phi = dist_phi(rng);
        scalar theta = dist_theta(rng);

        scalar sin_phi = std::sin(phi);
        scalar cos_phi = std::cos(phi);
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        point3 dir = vector3{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};
        auto hit_cocylindrical = cci.intersect(identity, ori, dir, cylindrical2, cylinder);
        if (hit_cocylindrical.status == intersection_status::e_inside)
        {
            if (write_files)
            {
                search_points
                    << hit_cocylindrical.point3[0] << ", "
                    << hit_cocylindrical.point3[1] << ", "
                    << hit_cocylindrical.point3[2] << "\n";
            }

            auto transform_indices = inner_finder(hit_cocylindrical.point2.value(), {1, 1});

            // Now loop over the local candidates: test case
            for (auto tfi : transform_indices)
            {
                auto hit_plane = pi.intersect(barrel_transforms[tfi], hit_cocylindrical.point3, dir, cart2, rect);
                if (hit_plane.status == intersection_status::e_inside)
                {
                    ++hit_modules_found;
                    if (write_files)
                    {
                        confirmed_points
                            << hit_plane.point3[0] << ", "
                            << hit_plane.point3[1] << ", "
                            << hit_plane.point3[2] << "\n";
                    }
                }
            }

            // Now loop over all candidates: heuristic case, cross-check
            if (test_against_all)
            {
                for (unsigned int itf = 0; itf < barrel_transforms.size(); ++itf)
                {
                    auto hit_plane_try = pi.intersect(barrel_transforms[itf], hit_cocylindrical.point3, dir, cart2, rect);
                    if (hit_plane_try.status == intersection_status::e_inside and hit_plane_try.path > 0.)
                    {
                        ++hit_modules_heuristic;
                        if (write_files)
                        {
                            heuristic_points
                                << hit_plane_try.point3[0] << ", "
                                << hit_plane_try.point3[1] << ", "
                                << hit_plane_try.point3[2] << "\n";
                        }
                    }
                }
            }
        }
    }

    if (write_files)
    {
        search_points.close();
        confirmed_points.close();
        heuristic_points.close();
    }

    if (test_against_all)
    {
        ASSERT_EQ(hit_modules_found, hit_modules_heuristic);
    }
}

scalar inner_r = 32.;
scalar outer_r = 44.;
scalar pos_z = 550.;
scalar stagger_z = 2.;
unsigned int n_phi_ec = 12;
scalar overlap_phi = 0.15;
scalar volume_inner_r = 29.;
scalar volume_outer_r = 50.;
scalar volume_min_z = 540.;
scalar volume_max_z = 560.;
unsigned int transform_offset = 0;

// Retrieving the endcap components for testing
auto endcap_components = detray::create_endcap_components(inner_r, outer_r, pos_z, stagger_z, n_phi_ec, overlap_phi,
                                                          volume_inner_r, volume_outer_r, volume_min_z, volume_max_z, transform_offset);

auto endcap_module = std::get<0>(endcap_components);
auto endcap_transforms = std::get<1>(endcap_components);
auto endcap_finders = std::get<2>(endcap_components);

// Test the basic barrel components function
TEST(__plugin, endcap_component_construction)
{
    // Rectangular module parameters
    ASSERT_EQ(endcap_module.size(), 3u);
    // Number of Barrel transforms
    ASSERT_EQ(endcap_transforms.size(), n_phi_ec);
    // Barrel finders
    ASSERT_EQ(endcap_finders.size(), 4u);

    // Check if the finders make sense
    auto inner_finder = endcap_finders[0];
    auto outer_finder = endcap_finders[1];
    auto nec_finder = endcap_finders[2];
    auto pec_finder = endcap_finders[3];

    using cylinder_grid = grid2<replace_populator<>, axis::circular<>, axis::closed<>, serializer2>;
    using disc_grid = grid2<replace_populator<>, axis::closed<>, axis::circular<>, serializer2>;

    using cylinder_finder = local_zone_finder<cylinder_grid>;
    using disc_finder = local_zone_finder<disc_grid>;

    auto inner_local_finder = inner_finder.target<cylinder_finder>();
    auto outer_local_finder = outer_finder.target<cylinder_finder>();
    auto nec_local_finder = nec_finder.target<disc_finder>();
    auto pec_local_finder = pec_finder.target<disc_finder>();

    EXPECT_TRUE(inner_local_finder != nullptr);
    EXPECT_TRUE(outer_local_finder != nullptr);
    EXPECT_TRUE(nec_local_finder != nullptr);
    EXPECT_TRUE(pec_local_finder != nullptr);

    const auto &inner_grid = inner_local_finder->grid();
    auto inner_grid_axis_p0 = inner_grid.axis_p0();
    auto inner_grid_axis_p1 = inner_grid.axis_p1();
    ASSERT_EQ(inner_grid_axis_p0.bins, n_phi_ec);
    ASSERT_EQ(inner_grid_axis_p1.bins, 1u);

    const auto &outer_grid = outer_local_finder->grid();
    auto outer_grid_axis_p0 = outer_grid.axis_p0();
    auto outer_grid_axis_p1 = outer_grid.axis_p1();
    ASSERT_EQ(outer_grid_axis_p0.bins, n_phi_ec);
    ASSERT_EQ(outer_grid_axis_p1.bins, 1u);

    const auto &ndisc_grid = nec_local_finder->grid();
    auto ndisc_grid_axis_p0 = ndisc_grid.axis_p0();
    auto ndisc_grid_axis_p1 = ndisc_grid.axis_p1();
    ASSERT_EQ(ndisc_grid_axis_p0.bins, 1u);
    ASSERT_EQ(ndisc_grid_axis_p1.bins, n_phi_ec);

    const auto &pdisc_grid = pec_local_finder->grid();
    auto pdisc_grid_axis_p0 = pdisc_grid.axis_p0();
    auto pdisc_grid_axis_p1 = pdisc_grid.axis_p1();
    ASSERT_EQ(pdisc_grid_axis_p0.bins, 1u);
    ASSERT_EQ(pdisc_grid_axis_p1.bins, n_phi_ec);

    using cyl_point2 = __plugin::cylindrical2::point2;
    using pol_point2 = __plugin::polar2::point2;

    scalar phi_step = 2 * M_PI / n_phi_ec;
    scalar r = 0.5 * (inner_r + outer_r);
    unsigned int imod = 0;

    // Test the entries of the inner/outer Barrel grids
    for (unsigned int iphi = 0; iphi < n_phi_ec; ++iphi)
    {
        scalar phi = -M_PI + iphi * phi_step;
        pol_point2 pol2{r, phi};
        cyl_point2 inner_p2(volume_inner_r * phi, pos_z);
        cyl_point2 outer_p2(volume_outer_r * phi, pos_z);

        auto inner_zone = inner_finder(inner_p2, {1, 1});
        auto outer_zone = outer_finder(outer_p2, {1, 1});
        auto nec_zone = nec_finder(pol2, {1, 1});
        auto pec_zone = pec_finder(pol2, {1, 1});

        ASSERT_EQ(inner_zone.size(), 3u);
        ASSERT_EQ(outer_zone.size(), 3u);
        ASSERT_EQ(nec_zone.size(), 3u);
        ASSERT_EQ(pec_zone.size(), 3u);

        ASSERT_EQ(inner_zone, outer_zone);
        ASSERT_EQ(nec_zone, pec_zone);
        ASSERT_EQ(inner_zone, nec_zone);
    }
}

TEST(__plugin, endcap_object_finding)
{

    std::ofstream search_points;
    std::ofstream heuristic_points;
    std::ofstream confirmed_points;

    if (write_files)
    {
        search_points.open("endcap_layer_search_points.csv");
        search_points << "x,y,z"
                      << "\n";

        confirmed_points.open("endcap_layer_confirmed_points.csv");
        confirmed_points << "x,y,z"
                         << "\n";

        heuristic_points.open("endcap_layer_heuristic_points.csv");
        heuristic_points << "x,y,z"
                         << "\n";
    }

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<scalar> dist_phi(-M_PI, M_PI);

    scalar min_theta = std::atan2(volume_inner_r, volume_min_z);
    scalar max_theta = std::atan2(volume_outer_r, volume_min_z);

    std::uniform_real_distribution<scalar> dist_theta(min_theta, max_theta);

    transform3 shift(vector3{0., 0., volume_min_z});
    ring2<> disc = {volume_inner_r, volume_outer_r};

    // Cylindrical and Cartesian local frames
    __plugin::cylindrical2 cylindrical2;
    __plugin::cartesian2 cart2;
    __plugin::polar2 pol2;

    vector3 ori = {0., 0., 0.};

    // Check if the finders make sense
    auto nec_finder = endcap_finders[2];

    trapezoid2<> trap = {endcap_module[0], endcap_module[1], endcap_module[2]};

    // Cross check
    unsigned int hit_modules_found = 0;
    unsigned int hit_modules_heuristic = 0;

    for (unsigned int itest = 0; itest < tests; ++itest)
    {

        scalar phi = dist_phi(rng);
        scalar theta = dist_theta(rng);

        scalar sin_phi = std::sin(phi);
        scalar cos_phi = std::cos(phi);
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        point3 dir = vector3{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};
        auto hit_disc = pi.intersect(shift, ori, dir, pol2, disc);

        if (hit_disc.status == intersection_status::e_inside)
        {
            if (write_files)
            {
                search_points
                    << hit_disc.point3[0] << ", "
                    << hit_disc.point3[1] << ", "
                    << hit_disc.point3[2] << "\n";
            }
        }

        auto transform_indices = nec_finder(hit_disc.point2.value(), {1, 1});

        // Now loop over the local candidates: test case
        for (auto tfi : transform_indices)
        {
            auto hit_plane = pi.intersect(endcap_transforms[tfi], hit_disc.point3, dir, cart2, trap);
            if (hit_plane.status == intersection_status::e_inside)
            {
                ++hit_modules_found;
                if (write_files)
                {
                    confirmed_points
                        << hit_plane.point3[0] << ", "
                        << hit_plane.point3[1] << ", "
                        << hit_plane.point3[2] << "\n";
                }
            }
        }

        // Now loop over all candidates: heuristic case, cross-check
        if (test_against_all)
        {
            for (unsigned int itf = 0; itf < endcap_transforms.size(); ++itf)
            {
                auto hit_plane_try = pi.intersect(endcap_transforms[itf], hit_disc.point3, dir, cart2, trap);
                if (hit_plane_try.status == intersection_status::e_inside and hit_plane_try.path > 0.)
                {
                    ++hit_modules_heuristic;
                    if (write_files)
                    {
                        heuristic_points
                            << hit_plane_try.point3[0] << ", "
                            << hit_plane_try.point3[1] << ", "
                            << hit_plane_try.point3[2] << "\n";
                    }
                }
            }
        }
    }

    if (write_files)
    {
        search_points.close();
        confirmed_points.close();
        heuristic_points.close();
    }

    if (test_against_all)
    {
        ASSERT_EQ(hit_modules_found, hit_modules_heuristic);
    }
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
