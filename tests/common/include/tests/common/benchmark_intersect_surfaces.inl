/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/cylinder3.hpp"
#include "masks/rectangle2.hpp"
#include "masks/single3.hpp"
#include "tests/common/test_surfaces.hpp"
#include "tools/concentric_cylinder_intersector.hpp"
#include "tools/cylinder_intersector.hpp"
#include "tools/planar_intersector.hpp"

#include <fstream>

#include <benchmark/benchmark.h>

using namespace detray;

using transform3 = __plugin::transform3;
using point3 = transform3::point3;
using vector3 = transform3::vector3;
using context = transform3::context;
using plane_surface = surface<transform3>;
__plugin::cartesian2 cart2;
__plugin::cylindrical2 cyl2;

unsigned int theta_steps = 1000;
unsigned int phi_steps = 1000;

bool screen_output = false;

dvector<scalar> dists = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

auto planes = planes_along_direction(dists, vector::normalize(vector3{1., 1., 1.}));

namespace __plugin
{
    // This test runs intersection with all surfaces of the TrackML detector
    static void BM_INTERSECT_PLANES(benchmark::State &state)
    {

        unsigned int sfhit = 0;
        unsigned int sfmiss = 0;

        planar_intersector pi;
        rectangle2<scalar> rect = {10., 20.};

        context ctx;

        point3 ori = {0., 0., 0.};

        for (auto _ : state)
        {
            // Loops of theta values
            for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
            {
                scalar theta = 0.85 + itheta * 0.2 / theta_steps;
                double sin_theta = std::sin(theta);
                double cos_theta = std::cos(theta);

                // Loops of phi values
                for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
                {
                    scalar phi = 0.7 + iphi * 0.2 / phi_steps;
                    double sin_phi = std::sin(phi);
                    double cos_phi = std::cos(phi);

                    vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

                    for (auto plane : planes)
                    {
                        auto is = pi.intersect(plane, ori, dir, ctx, cart2, rect);
                        if (is._status == e_inside)
                        {
                            ++sfhit;
                        }
                        else
                        {
                            ++sfmiss;
                        }
                    }
                }
            }
        }

        if (screen_output)
        {
            std::cout << "Surfaces inside / outside : " << sfhit << " / " << sfmiss << std::endl;
        }
    }

    BENCHMARK(BM_INTERSECT_PLANES);

    // This test runs intersection with all surfaces of the TrackML detector
    static void BM_INTERSECT_CYLINDERS(benchmark::State &state)
    {

        unsigned int sfhit = 0;
        unsigned int sfmiss = 0;

        cylinder_intersector ci;
        dvector<cylinder3<scalar>> cylinders;
        for (auto r : dists)
        {
            cylinders.push_back(cylinder3<scalar>{r, 10.});
        }

        context ctx;
        surface<transform3> plain(std::move(transform3()), 0, 0);

        point3 ori = {0., 0., 0.};

        for (auto _ : state)
        {
            // Loops of theta values
            for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
            {
                scalar theta = 0.85 + itheta * 0.2 / theta_steps;
                double sin_theta = std::sin(theta);
                double cos_theta = std::cos(theta);

                // Loops of phi values
                for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
                {
                    scalar phi = 0.7 + iphi * 0.2 / phi_steps;
                    double sin_phi = std::sin(phi);
                    double cos_phi = std::cos(phi);

                    vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

                    for (auto cylinder : cylinders)
                    {
                        auto is = ci.intersect(plain, ori, dir, ctx, cyl2, cylinder);
                        if (is._status == e_inside)
                        {
                            ++sfhit;
                        }
                        else
                        {
                            ++sfmiss;
                        }
                    }
                }
            }
        }

        if (screen_output)
        {
            std::cout << "Surfaces inside / outside : " << sfhit << " / " << sfmiss << std::endl;
        }
    }

    BENCHMARK(BM_INTERSECT_CYLINDERS);

    // This test runs intersection with all surfaces of the TrackML detector
    static void BM_INTERSECT_CONCETRIC_CYLINDERS(benchmark::State &state)
    {
        unsigned int sfhit = 0;
        unsigned int sfmiss = 0;

        concentric_cylinder_intersector cci;
        dvector<single3<scalar, 2>> cylinder_halfs;
        for (auto r : dists)
        {
            cylinder_halfs.push_back(single3<scalar, 2>{10.});
        }

        context ctx;
        surface<transform3> plain(std::move(transform3()), 0, 0);

        point3 ori = {0., 0., 0.};

        for (auto _ : state)
        {
            // Loops of theta values
            for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
            {
                scalar theta = 0.85 + itheta * 0.2 / theta_steps;
                double sin_theta = std::sin(theta);
                double cos_theta = std::cos(theta);

                // Loops of phi values
                for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
                {
                    scalar phi = 0.7 + iphi * 0.2 / phi_steps;
                    double sin_phi = std::sin(phi);
                    double cos_phi = std::cos(phi);

                    vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

                    for (auto [ic, cylinder_half] : enumerate(cylinder_halfs))
                    {
                        auto is = cci.intersect(plain, dists[ic], ori, dir, ctx, cyl2, cylinder_half);
                        if (is._status == e_inside)
                        {
                            ++sfhit;
                        }
                        else
                        {
                            ++sfmiss;
                        }
                    }
                }
            }
        }

        if (screen_output)
        {
            std::cout << "Surfaces inside / outside : " << sfhit << " / " << sfmiss << std::endl;
        }
    }

    BENCHMARK(BM_INTERSECT_CONCETRIC_CYLINDERS);

} // namespace __plugin

BENCHMARK_MAIN();
