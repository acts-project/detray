/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "tools/intersection_kernel.hpp"
#include "tools/line_stepper.hpp"
#include "tools/navigator.hpp"
#include "tools/propagator.hpp"

#include <iostream>
#include <fstream>

#include <benchmark/benchmark.h>

using namespace detray;

using transform3 = __plugin::transform3;
using point3 = point3;
using vector3 = vector3;
using surface = surface_base<transform3>;

__plugin::cartesian2 cartesian2;
using point2 = __plugin::point2;

unsigned int theta_steps = 100;
unsigned int phi_steps = 100;
bool stream_file = false;

/** Read the detector from file */
auto read_detector()
{
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr)
    {
        throw std::ios_base::failure("Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string name = "tml";
    std::string surfaces = data_directory + "tml.csv";
    std::string grids = data_directory + "tml-surface-grids.csv";
    std::string volumes = data_directory + "tml-layer-volumes.csv";
    return detray::detector_from_csv<static_transform_store>(name, surfaces, grids, volumes);
};

auto d = read_detector();

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

        unsigned int hits = 0;
        unsigned int missed = 0;

        point3 ori = {0., 0., 0.};

        for (auto _ : state)
        {

            track<static_transform_store::context> track;
            track.pos = point3{0., 0., 0.};

            // Loops of theta values
            for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
            {
                scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
                scalar sin_theta = std::sin(theta);
                scalar cos_theta = std::cos(theta);

                // Loops of phi values
                for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
                {
                    // The direction
                    scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
                    scalar sin_phi = std::sin(phi);
                    scalar cos_phi = std::cos(phi);
                    track.dir = {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

                    // Loop over volumes
                    for (const auto &v : d.volumes())
                    {
                        const auto &surfaces = v.surfaces();

                        // Loop over surfaces
                        for (const auto &s : surfaces.objects())
                        {
                            auto sfi_surface = intersect(track, s, surfaces.transforms(), surfaces.masks());

                            const auto &sfi = std::get<0>(sfi_surface);
                            if (sfi.status == intersection_status::e_inside)
                            {
                                if (stream_file)
                                {
                                    hit_out << sfi.p3[0] << "," << sfi.p3[1] << "," << sfi.p3[2] << "\n";
                                }
                                ++hits;
                            }
                            else
                            {
                                ++missed;
                            }
                        }
                    }
                }
            }
        }

        std::cout << "[detray] hits / missed / total = " << hits << " / " << missed << " / " << hits + missed << std::endl;
        if (stream_file)
        {
            hit_out.close();
        }
    }

    BENCHMARK(BM_INTERSECT_ALL);

    // This test runs intersection with all surfaces using the propagator
    static void BM_PROPAGATE_ALL(benchmark::State &state)
    {

        unsigned int hits = 0;
        unsigned int missed = 0;

        point3 ori = {0., 0., 0.};

        // Create the navigator
        using detray_navigator = navigator<decltype(d)>;
        using detray_track = track<static_transform_store::context>;

        detray_track traj;
        traj.pos = {0., 0., 0.};
        traj.ctx = static_transform_store::context{};
        traj.momentum = 100.;
        traj.overstep_tolerance = -1e-4;
        traj.path_limit = 1000.;

        using detray_stepper = line_stepper<detray_track>;

        detray_stepper s;
        detray_navigator n(std::move(d));

        using detray_propagator = propagator<detray_stepper, detray_navigator>;
        detray_propagator p(std::move(s), std::move(n));

        void_track_inspector vi;

        for (auto _ : state)
        {

            // Loops of theta values
            for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
            {
                scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
                scalar sin_theta = std::sin(theta);
                scalar cos_theta = std::cos(theta);

                // Loops of phi values
                for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
                {
                    // The direction
                    scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
                    scalar sin_phi = std::sin(phi);
                    scalar cos_phi = std::cos(phi);
                    traj.dir = {1., 0., 0. }; // {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

                    auto end = p.propagate(traj, vi);


                }
            }
        }
    }

    BENCHMARK(BM_PROPAGATE_ALL);

} // namespace __plugin

BENCHMARK_MAIN();
