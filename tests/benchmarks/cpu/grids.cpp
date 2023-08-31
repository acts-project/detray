/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray core include(s).
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"

// Detray test include(s).
#include "detray/test/types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Google benchmark include(s).
#include <benchmark/benchmark.h>

// Use the detray:: namespace implicitly.
using namespace detray;

/*namespace {

/// Make a regular grid for the tests.
auto make_regular_grid(vecmem::memory_resource &mr) {

    // Return the grid.
    return grid2<replace_populator, axis::regular, axis::regular, serializer2>{
        {25u, 0.f, 25.f, mr}, {60u, 0.f, 60.f, mr}, mr};
}

}  // namespace

// This runs a reference test with a regular grid structure
void BM_REGULAR_GRID_BIN(benchmark::State &state) {

    // Set up the tested grid object.
    vecmem::host_memory_resource host_mr;
    auto g2r = make_regular_grid(host_mr);

    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2 p = {static_cast<scalar>((rand() % 50)) * 0.5f,
                              static_cast<scalar>((rand() % 120)) * 0.5f};
            g2r.bin(p);
        }
    }
}

BENCHMARK(BM_REGULAR_GRID_BIN)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

void BM_REGULAR_GRID_ZONE(benchmark::State &state) {

    // Set up the tested grid object.
    vecmem::host_memory_resource host_mr;
    auto g2r = make_regular_grid(host_mr);

    // Helper zone object.
    static const darray<dindex, 2> zone22 = {2u, 2u};

    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2 p = {static_cast<scalar>((rand() % 50)) * 0.5f,
                              static_cast<scalar>((rand() % 120)) * 0.5f};
            g2r.zone(p, {zone22, zone22});
        }
    }
}

BENCHMARK(BM_REGULAR_GRID_ZONE)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

namespace {

/// Make an irregular grid for the tests.
auto make_irregular_grid(vecmem::memory_resource &mr) {

    // Fill a 25 x 60 grid with an "irregular" axis
    dvector<scalar> xboundaries, yboundaries;
    xboundaries.reserve(25u);
    yboundaries.reserve(60u);
    for (scalar i = 0.f; i < 61.f; i += 1.f) {
        if (i < 26.f) {
            xboundaries.push_back(i);
        }
        yboundaries.push_back(i);
    }

    return grid2<replace_populator, axis::irregular, axis::irregular,
                 serializer2>({xboundaries, mr}, {yboundaries, mr}, mr);
}

}  // namespace

// This runs a reference test with a irregular grid structure
void BM_IRREGULAR_GRID_BIN(benchmark::State &state) {

    // Set up the tested grid object.
    vecmem::host_memory_resource host_mr;
    auto g2irr = make_irregular_grid(host_mr);

    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2 p = {static_cast<scalar>((rand() % 50)) * 0.5f,
                              static_cast<scalar>((rand() % 120)) * 0.5f};
            g2irr.bin(p);
        }
    }
}

BENCHMARK(BM_IRREGULAR_GRID_BIN)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

void BM_IRREGULAR_GRID_ZONE(benchmark::State &state) {

    // Set up the tested grid object.
    vecmem::host_memory_resource host_mr;
    auto g2irr = make_irregular_grid(host_mr);

    // Helper zone object.
    static const darray<dindex, 2> zone22 = {2u, 2u};

    for (auto _ : state) {
        for (unsigned int itest = 0u; itest < 1000000u; ++itest) {
            test::point2 p = {static_cast<scalar>((rand() % 50)) * 0.5f,
                              static_cast<scalar>((rand() % 120)) * 0.5f};
            g2irr.zone(p, {zone22, zone22});
        }
    }
}

BENCHMARK(BM_IRREGULAR_GRID_ZONE)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
*/
