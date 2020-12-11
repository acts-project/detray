/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/intersection.hpp"
#include "tools/planar_intersector.hpp"
#include "tests/io/read_csv.hpp"

#include <fstream>

#include <benchmark/benchmark.h>

using namespace detray;

using transform3 = __plugin::transform3;
using point3 = transform3::point3;
using vector3 = transform3::vector3;
using context = transform3::context;
using surface = surface<transform3>;

__plugin::cartesian2 cartesian2;
using point2 = __plugin::cartesian2::point2;

using planar_intersection = intersection<scalar,point3,point2>;

unsigned int theta_steps = 100;
unsigned int phi_steps = 100;
bool stream_file = false;


/** Read the detector from file */
auto read_from_file()
{
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr)
    {
        throw std::ios_base::failure("Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);
    return read_csv<transform3>(data_directory + std::string("/tml-detector.csv"));
};

auto detector = read_from_file();

namespace __plugin
{
    // This test runs intersection with all surfaces of the TrackML detector
    static void BM_INTERSECT_ALL(benchmark::State &state)
    {

        std::ofstream hit_out;
        if (stream_file)
        {
            hit_out.open("tml_hits.csv");
        }

        unsigned int hits_inside = 0;
        unsigned int hits_outside = 0;
        unsigned int hits_missed = 0;
        context ctx;

        point3 ori = {0., 0., 0.};

        for (auto _ : state)
        {
            // Loops of theta values
            for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
            {
                scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
                double sin_theta = std::sin(theta);
                double cos_theta = std::cos(theta);

                // Loops of phi values
                for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
                {
                    scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
                    double sin_phi = std::sin(phi);
                    double cos_phi = std::cos(phi);

                    vector3 dir = {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

                    for (auto &volume : detector.volumes)
                    {
                        for (auto &layer : volume.layers)
                        {
                            for (auto &surface : layer.surfaces)
                            {
                             
                                planar_intersection hit;

                                auto group_index = surface.mask()[0];                                
                                auto mask_index = surface.mask()[1];
                                if (group_index == 0){
                                    auto mask = std::get<0>(layer.masks)[mask_index];
                                    hit = mask.intersector().intersect(surface, ori, dir, ctx, cartesian2, mask);
                                } else {
                                    auto mask = std::get<1>(layer.masks)[mask_index];
                                    hit = mask.intersector().intersect(surface, ori, dir, ctx, cartesian2, mask);
                                }

                                if (hit._status == e_inside)
                                {
                                    ++hits_inside;
                                    if (stream_file)
                                    {
                                        hit_out << hit._point3[0] << "," << hit._point3[1] << "," << hit._point3[2] << "\n";
                                    }
                                }
                                else if (hit._status == e_outside)
                                {
                                    ++hits_outside;
                                }
                                else
                                {
                                    ++hits_missed;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (stream_file)
        {
            hit_out.close();
        }
    }

    BENCHMARK(BM_INTERSECT_ALL);

} // namespace __plugin

BENCHMARK_MAIN();
