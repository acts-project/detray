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

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "tests/common/read_geometry.hpp"

using namespace detray;

#ifdef DETRAY_BENCHMARKS_REP
unsigned int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
unsigned int gbench_repetitions = 0;
#endif

auto [d, name_map] = detray_tests::read_from_csv(detray_tests::tml_files);
const unsigned int itest = 10000;

namespace __plugin {
// This test a reference run to deduce the random number
static void BM_FIND_VOLUMES(benchmark::State &state) {
    auto volume_grid = d.volume_search_grid();

    const auto &axis0 = volume_grid.axis_p0();
    const auto &axis1 = volume_grid.axis_p1();

    auto range0 = axis0.span();
    auto range1 = axis1.span();

    scalar step0 = (range0[1] - range0[0]) / itest;
    scalar step1 = (range0[1] - range0[0]) / itest;

    size_t successful = 0;
    size_t unsuccessful = 0;

    for (auto _ : state) {
        for (unsigned int i1 = 0; i1 < itest; ++i1) {
            for (unsigned int i0 = 0; i0 < itest; ++i0) {
                vector3 rz{i0 * step0, 0., i1 * step1};
                auto &v = d.volume_by_pos(rz);

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

}  // namespace __plugin

BENCHMARK_MAIN();
