/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <vecmem/memory/host_memory_resource.hpp>
// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "detray/definitions/indexing.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"

using namespace detray;
using namespace __plugin;

namespace {
// memory resource
vecmem::host_memory_resource host_mr;

darray<dindex, 2> zone22 = {2u, 2u};

// TrackML detector has 25 x 60 cells int he detector grid
using grid2r =
    grid2<replace_populator, axis::regular, axis::regular, serializer2>;
grid2r::axis_p0_type xaxisr{25u, 0.f, 25.f, host_mr};
grid2r::axis_p1_type yaxisr{60u, 0.f, 60.f, host_mr};

grid2r g2r(std::move(xaxisr), std::move(yaxisr), host_mr);

// This runs a reference test with a regular grid structure
static void BM_REGULAR_GRID_BIN(benchmark::State &state) {
    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2<detray::scalar> p = {
                static_cast<scalar>((rand() % 50)) * 0.5f,
                static_cast<scalar>((rand() % 120)) * 0.5f};
            g2r.bin(p);
        }
    }
}

static void BM_REGULAR_GRID_ZONE(benchmark::State &state) {
    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2<detray::scalar> p = {
                static_cast<scalar>((rand() % 50)) * 0.5f,
                static_cast<scalar>((rand() % 120)) * 0.5f};
            g2r.zone(p, {zone22, zone22});
        }
    }
}

auto construct_irregular_grid() {
    // Fill a 25 x 60 grid with an "irregular" axis
    dvector<scalar> xboundaries = {};
    xboundaries.reserve(25u);
    dvector<scalar> yboundaries = {};
    yboundaries.reserve(60u);

    for (unsigned int i = 0u; i < 61u; ++i) {
        if (i < 26u) {
            xboundaries.push_back(i);
        }
        yboundaries.push_back(i);
    }

    using grid2ir =
        grid2<replace_populator, axis::irregular, axis::irregular, serializer2>;

    grid2ir::axis_p0_type xaxisir{xboundaries, host_mr};
    grid2ir::axis_p1_type yaxisir{yboundaries, host_mr};

    return grid2ir(std::move(xaxisir), std::move(yaxisir), host_mr);
}

auto g2irr = construct_irregular_grid();

// This runs a reference test with a irregular grid structure
static void BM_IRREGULAR_GRID_BIN(benchmark::State &state) {
    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2<detray::scalar> p = {
                static_cast<scalar>((rand() % 50)) * 0.5f,
                static_cast<scalar>((rand() % 120)) * 0.5f};
            g2irr.bin(p);
        }
    }
}

static void BM_IRREGULAR_GRID_ZONE(benchmark::State &state) {
    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2<detray::scalar> p = {
                static_cast<scalar>((rand() % 50)) * 0.5f,
                static_cast<scalar>((rand() % 120)) * 0.5f};
            g2irr.zone(p, {zone22, zone22});
        }
    }
}

// BENCHMARK(BM_RERERENCE_GRID);
BENCHMARK(BM_REGULAR_GRID_BIN);
BENCHMARK(BM_REGULAR_GRID_ZONE);
BENCHMARK(BM_IRREGULAR_GRID_BIN);
BENCHMARK(BM_IRREGULAR_GRID_ZONE);

}  // namespace

BENCHMARK_MAIN();
