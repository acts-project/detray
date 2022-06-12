/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/detector.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/propagator/track.hpp"
#include "detray/utils/enumerate.hpp"
#include "tests/common/tools/detector_metadata.hpp"
#include "tests/common/tools/read_geometry.hpp"
#include "tests/common/tools/track_generators.hpp"

using namespace detray;

#ifdef DETRAY_BENCHMARKS_REP
unsigned int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
unsigned int gbench_repetitions = 0;
#endif

unsigned int theta_steps = 100;
unsigned int phi_steps = 100;
bool stream_file = false;

vecmem::host_memory_resource host_mr;
auto [d, name_map] =
    read_from_csv<detector_registry::tml_detector>(tml_files, host_mr);

using detector_t = decltype(d);
constexpr auto k_surfaces = detector_t::objects::e_surface;

using detray_context = detector_t::context;
detray_context default_context;

const auto data_core = d.data(default_context);

namespace __plugin {
// This test runs intersection with all surfaces of the TrackML detector
static void BM_INTERSECT_ALL(benchmark::State &state) {

    /*std::ofstream hit_out;
    if (stream_file)
    {
        hit_out.open("tml_hits.csv");
    }*/
    unsigned int hits = 0;
    unsigned int missed = 0;

    // point3 ori = {0., 0., 0.};

    for (auto _ : state) {
        point3<detray::scalar> pos{0., 0., 0.};

        // Iterate through uniformly distributed momentum directions
        for (const auto track : uniform_track_generator<free_track_parameters>(
                 theta_steps, phi_steps, pos)) {

            // Loop over volumes
            for (const auto &v : d.volumes()) {
                // Loop over all surfaces in volume
                for (const auto sf : range(data_core.surfaces, v)) {
                    auto sfi = intersect(track, sf, data_core.transforms,
                                         data_core.masks);

                    benchmark::DoNotOptimize(hits);
                    benchmark::DoNotOptimize(missed);
                    if (sfi.status == intersection::status::e_inside) {
                        /* state.PauseTiming();
                        if (stream_file)
                        {
                            hit_out << sfi.p3[0] << "," << sfi.p3[1] << ","
                        << sfi.p3[2]
                        << "\n";
                        }
                        state.ResumeTiming();*/
                        ++hits;
                    } else {
                        ++missed;
                    }
                }
            }
        }
    }

#ifndef DETRAY_BENCHMARKS_MULTITHREAD
    std::cout << "[detray] hits / missed / total = " << hits << " / " << missed
              << " / " << hits + missed << std::endl;
#endif
    /**if (stream_file)
    {
        hit_out.close();
    }*/
}

BENCHMARK(BM_INTERSECT_ALL)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

}  // namespace __plugin

BENCHMARK_MAIN();
