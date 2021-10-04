/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "tools/intersection_kernel.hpp"

using namespace detray;

/** Read the detector from file */
auto read_detector() {
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure(
            "Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string name = "tml";
    std::string surfaces = data_directory + "tml.csv";
    std::string volumes = data_directory + "tml-layer-volumes.csv";
    std::string grids = data_directory + "tml-surface-grids.csv";
    std::string grid_entries = "";
    //std::map<dindex, std::string> name_map{};

    /*return detray::detector_from_csv<>(name, surfaces, volumes, grids,
                                       grid_entries, name_map);*/
    return detray::detector_from_csv<>(name, surfaces, volumes, grids,
                                       grid_entries);
};

/** Return a random ray within some range */
auto shoot_ray(std::pair<point3, point3> origin,
               std::pair<point3, point3> direction) {
    track<static_transform_store<>::context> ray;
    ray.pos = origin.first;
    ray.dir = direction.first;

    return ray;
};

unsigned int theta_steps = 100;
unsigned int phi_steps = 100;
const unsigned int itest = 10000;

auto d = read_detector();
const auto &portals = d.portals();
const auto &masks = d.masks();
using links_type = typename decltype(d)::geometry::portal_links;

namespace __plugin {
// This test runs intersection with all surfaces of the TrackML detector
TEST(ALGEBRA_PLUGIN, random_rays) {
    /*std::ofstream hit_out;
    if (stream_file)
    {
        hit_out.open("tml_hits.csv");
    }*/
    unsigned int hits = 0;
    unsigned int missed = 0;

    using detray_context = decltype(d)::transform_store::context;
    detray_context default_context;

    const auto ori = std::make_pair<point3, point3>({0., 0., 0.}, {0., 0., 0.});

    // Loops of theta values
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            const auto dir = std::make_pair<point3, point3>(
                {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta},
                {0., 0., 0.});

            const auto ray = shoot_ray(ori, dir);
            std::vector<std::pair<dindex, intersection>> volume_record;
            // Loop over volumes
            for (const auto &v : d.volumes()) {
                // Loop over portals
                const auto &pt_range = v.template range<decltype(d)::objects::e_portal>();
                links_type links = {};
                for (size_t pi = pt_range[0]; pi < pt_range[1]; pi++) {
                    auto sfi = intersect(
                        ray, portals[pi],
                        d.transforms(default_context),
                        masks, links);

                    if (sfi.status == intersection_status::e_inside) {
                        /* state.PauseTiming();
                        if (stream_file)
                        {
                            hit_out << sfi.p3[0] << "," << sfi.p3[1] << "," <<
                        sfi.p3[2] << "\n";
                        }
                        state.ResumeTiming();*/
                        volume_record.emplace_back(v.index(), sfi);
                    }
                }
            }
            std::cout << "ray o (" << ori.first[0] << ", " << ori.first[1]
                      << ", " << ori.first[2] << ") -> (" << dir.first[0]
                      << ", " << dir.first[1] << ", " << dir.first[2]
                      << "):" << std::endl;
            for (const auto &v_id : volume_record) {
                std::cout << v_id.first << " (" << v_id.second.path << ")"
                          << std::endl;
            }
        }
    }

    /**if (stream_file)
    {
        hit_out.close();
    }*/
}

// This test a reference run to deduce the random number
/*static void BM_G(benchmark::State &state)
{
    auto volume_grid = d.volume_search_grid();

    const auto &axis0 = volume_grid.axis_p0();
    const auto &axis1 = volume_grid.axis_p1();

    auto range0 = axis0.span();
    auto range1 = axis1.span();

    scalar step0 = (range0[1] - range0[0]) / itest;
    scalar step1 = (range0[1] - range0[0]) / itest;

    size_t successful = 0;
    size_t unsuccessful = 0;

    for (auto _ : state)
    {
        for (unsigned int i1 = 0; i1 < itest; ++i1)
        {
            for (unsigned int i0 = 0; i0 < itest; ++i0)
            {
                vector3 rz {i0 * step0, 0., i1 * step1};
                auto &v = d.indexed_volume(rz);

                benchmark::DoNotOptimize(successful);
                benchmark::DoNotOptimize(unsuccessful);
                if (v.index() == dindex_invalid)
                {
                    ++unsuccessful;
                }
                else
                {
                    ++successful;
                }
                benchmark::ClobberMemory();
            }
        }
    }


    #ifndef DETRAY_BENCHMARKS_MULTITHREAD
    std::cout << "Successful   : " << successful << std::endl;
    std::cout << "Unsuccessful : " << unsuccessful << std::endl;
    #endif
}*/

/*BENCHMARK(BM_FIND_VOLUMES)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
->ThreadPerCpu()
#endif
->Unit(benchmark::kMillisecond)
->Repetitions(gbench_repetitions)
->DisplayAggregatesOnly(true);*/

}  // namespace __plugin

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
