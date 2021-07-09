/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"

#include <benchmark/benchmark.h>

#include <string>
#include <iostream>

using namespace detray;

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
    std::string surfaces = data_directory + "odd.csv";
    std::string volumes = data_directory + "odd-layer-volumes.csv";
    std::string grids = data_directory + "odd-surface-grids.csv";
    std::string grid_entries = "";
    return detray::detector_from_csv<>(name, surfaces, volumes, grids, grid_entries);
};

auto d = read_detector();
const unsigned int itest = 10000;

namespace __plugin
{
    // This test a reference run to deduce the random number
    static void BM_FIND_VOLUMES(benchmark::State &state)
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

        for (unsigned int i1 = 0; i1 < itest; ++i1)
        {
            for (unsigned int i0 = 0; i0 < itest; ++i0)
            {
                vector3 rz {i0 * step0, 0., i1 * step1};
                auto &v = d.indexed_volume(rz);
                if (v.index() == dindex_invalid)
                {
                    ++unsuccessful;
                }
                else
                {
                    ++successful;
                }
            }
        }

        std::cout << "Successful   : " << successful << std::endl;
        std::cout << "Unsuccessful : " << unsuccessful << std::endl;
    }

    BENCHMARK(BM_FIND_VOLUMES);

} // namespace __plugin

BENCHMARK_MAIN();
