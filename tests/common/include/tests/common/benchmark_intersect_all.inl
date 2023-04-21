/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// System include(s)
#include <fstream>
#include <iostream>
#include <map>
#include <string>

using namespace detray;

#ifdef DETRAY_BENCHMARKS_REP
int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
int gbench_repetitions = 0;
#endif

unsigned int theta_steps{100u};
unsigned int phi_steps{100u};
bool stream_file = false;

// Detector configuration
constexpr std::size_t n_brl_layers{4u};
constexpr std::size_t n_edc_layers{7u};
vecmem::host_memory_resource host_mr;
auto d = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

using detector_t = decltype(d);
using intersection_t = intersection2D<typename detector_t::surface_type,
                                      typename detector_t::transform3>;

using detray_context = detector_t::geometry_context;
detray_context geo_context;

const auto &masks = d.mask_store();
const auto &transforms = d.transform_store(geo_context);

namespace __plugin {

// This test runs intersection with all surfaces of the TrackML detector
static void BM_INTERSECT_ALL(benchmark::State &state) {

    /*std::ofstream hit_out;
    if (stream_file)
    {
        hit_out.open("tml_hits.csv");
    }*/
    std::size_t hits{0u};
    std::size_t missed{0u};

    // point3 ori = {0.f, 0.f, 0.f};

    for (auto _ : state) {
        point3<detray::scalar> pos{0.f, 0.f, 0.f};
        std::vector<intersection_t> intersections{};
        std::size_t n_surfaces{0u};

        // Iterate through uniformly distributed momentum directions
        for (const auto track :
             uniform_track_generator<free_track_parameters<transform3<scalar>>>(
                 theta_steps, phi_steps, pos)) {

            // Loop over volumes
            for (const auto &v : d.volumes()) {
                // Loop over all surfaces in volume
                for (const auto &sf : d.surfaces(v)) {

                    masks.template visit<intersection_initialize>(
                        sf.mask(), intersections, detail::ray(track), sf,
                        transforms);

                    ++n_surfaces;

                    /* state.PauseTiming();
                    if (stream_file) {
                        for (const auto& sfi : intersections) {
                            hit_out << sfi.p3[0] << "," << sfi.p3[1] << ","
                                                 << sfi.p3[2] << "\n";
                        }
                    }
                    state.ResumeTiming();*/
                }
                benchmark::DoNotOptimize(hits);
                benchmark::DoNotOptimize(missed);

                hits += intersections.size();
                missed += n_surfaces - intersections.size();

                n_surfaces = 0u;
                intersections.clear();
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
