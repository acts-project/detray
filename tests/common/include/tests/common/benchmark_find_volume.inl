/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <iostream>
#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/detector_metadata.hpp"

using namespace detray;

namespace {

using transform3 = __plugin::transform3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;

}  // namespace

#ifdef DETRAY_BENCHMARKS_REP
unsigned int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
unsigned int gbench_repetitions = 0;
#endif

// Detector configuration
constexpr std::size_t n_brl_layers{4};
constexpr std::size_t n_edc_layers{7};
vecmem::host_memory_resource host_mr;
auto d = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

const unsigned int itest = 10000;

// Benchmarks the cost of searching a volume by position
static void BM_FIND_VOLUMES(benchmark::State &state) {
    auto &volume_grid = d.volume_search_grid();

    const auto &axis_r = volume_grid.get_axis<n_axis::label::e_r>();
    const auto &axis_z = volume_grid.get_axis<n_axis::label::e_z>();

    // Get a rough step size from irregular axes
    auto range0 = axis_r.span();
    auto range1 = axis_z.span();

    scalar step0{(range0[1] - range0[0]) / itest};
    scalar step1{(range1[1] - range1[0]) / itest};

    std::size_t successful{0};
    std::size_t unsuccessful{0};

    for (auto _ : state) {
        for (unsigned int i1 = 0; i1 < itest; ++i1) {
            for (unsigned int i0 = 0; i0 < itest; ++i0) {
                vector3 rz{i0 * step0, 0., i1 * step1};
                const auto &v = d.volume_by_pos(rz);

                benchmark::DoNotOptimize(successful);
                benchmark::DoNotOptimize(unsuccessful);
                if (v.index() == dindex_invalid) {
                    ++unsuccessful;
                } else {
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
}

BENCHMARK(BM_FIND_VOLUMES)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

BENCHMARK_MAIN();
