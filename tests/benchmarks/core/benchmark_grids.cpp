/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include <vecmem/memory/host_memory_resource.hpp>
#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "utils/indexing.hpp"

#include <benchmark/benchmark.h>

using namespace detray;
using namespace __plugin;

namespace
{
    // memory resource
    vecmem::host_memory_resource host_mr;
    
    darray<dindex, 2> zone22 = { 2u, 2u };

    // This runs a reference test with random numbers only
    static void BM_RERERENCE_GRID(benchmark::State &state)
    {
        for (auto _ : state)
        {
            for (unsigned int itest = 0; itest < 1000000; ++itest)
            {
                test::point2 p = {static_cast<scalar>((rand() % 50) * 0.5), static_cast<scalar>((rand() % 120) * 0.5)};
            }
        }
    }

    replace_populator<> replacer;
    serializer2 serializer;

    // TrackML detector has 25 x 60 cells int he detector grid
    axis::regular<> xaxisr = axis::regular<>{25, 0., 25.};
    axis::regular<> yaxisr = axis::regular<>{60, 0., 60.};
    using grid2r = grid2<decltype(replacer), decltype(xaxisr), decltype(yaxisr), decltype(serializer)>;
    grid2r g2r(std::move(xaxisr), std::move(yaxisr), host_mr);

    // This runs a reference test with a regular grid structure
    static void BM_REGULAR_GRID_BIN(benchmark::State &state)
    {
        for (auto _ : state)
        {
            for (unsigned int itest = 0; itest < 1000000; ++itest)
            {
                test::point2 p = {static_cast<scalar>((rand() % 50) * 0.5), static_cast<scalar>((rand() % 120) * 0.5)};
                g2r.bin(p);
            }
        }
    }

    static void BM_REGULAR_GRID_ZONE(benchmark::State &state)
    {
        for (auto _ : state)
        {
            for (unsigned int itest = 0; itest < 1000000; ++itest)
            {
                test::point2 p = {static_cast<scalar>((rand() % 50) * 0.5), static_cast<scalar>((rand() % 120) * 0.5)};
                g2r.zone(p, {zone22, zone22});
            }
        }
    }

    auto construct_irregular_grid()
    {
        // Fill a 25 x 60 grid with an "irregular" axis
        dvector<scalar> xboundaries = {};
        xboundaries.reserve(25);
        dvector<scalar> yboundaries = {};
        yboundaries.reserve(60);

        for (unsigned int i = 0; i < 61; ++i)
        {
            if (i < 26)
            {
                xboundaries.push_back(i);
            }
            yboundaries.push_back(i);
        }

        axis::irregular<> xaxisir{xboundaries};
        axis::irregular<> yaxisir{yboundaries};

        using grid2ir = grid2<decltype(replacer), decltype(xaxisir), decltype(yaxisir), decltype(serializer)>;

        return grid2ir(std::move(xaxisir), std::move(yaxisir), host_mr);
    }

    auto g2irr = construct_irregular_grid();

    // This runs a reference test with a irregular grid structure
    static void BM_IRREGULAR_GRID_BIN(benchmark::State &state)
    {
        for (auto _ : state)
        {
            for (unsigned int itest = 0; itest < 1000000; ++itest)
            {
                test::point2 p = {static_cast<scalar>((rand() % 50) * 0.5), static_cast<scalar>((rand() % 120) * 0.5)};
                g2irr.bin(p);
            }
        }
    }

    static void BM_IRREGULAR_GRID_ZONE(benchmark::State &state)
    {
        for (auto _ : state)
        {
            for (unsigned int itest = 0; itest < 1000000; ++itest)
            {
                test::point2 p = {static_cast<scalar>((rand() % 50) * 0.5), static_cast<scalar>((rand() % 120) * 0.5)};
                g2irr.zone(p, {zone22, zone22});
            }
        }
    }

    BENCHMARK(BM_RERERENCE_GRID);
    BENCHMARK(BM_REGULAR_GRID_BIN);
    BENCHMARK(BM_REGULAR_GRID_ZONE);
    BENCHMARK(BM_IRREGULAR_GRID_BIN);
    BENCHMARK(BM_IRREGULAR_GRID_ZONE);

} // namespace

BENCHMARK_MAIN();
