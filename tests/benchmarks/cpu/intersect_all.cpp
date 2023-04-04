/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// System include(s)
#include <iostream>
#include <map>
#include <string>

// Use the detray:: namespace implicitly.
using namespace detray;

// This test runs intersection with all surfaces of the TrackML detector
void BM_INTERSECT_ALL(benchmark::State &state) {

    static const unsigned int theta_steps{100u};
    static const unsigned int phi_steps{100u};

    // Detector configuration
    static constexpr std::size_t n_brl_layers{4u};
    static constexpr std::size_t n_edc_layers{7u};
    vecmem::host_memory_resource host_mr;
    auto d = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    decltype(d)::geometry_context geo_context;

    const auto &masks = d.mask_store();
    const auto &transforms = d.transform_store(geo_context);

    unsigned int hits{0u};
    unsigned int missed{0u};

    for (auto _ : state) {
        test::point3 pos{0.f, 0.f, 0.f};

        // Iterate through uniformly distributed momentum directions
        for (const auto track :
             uniform_track_generator<free_track_parameters<test::transform3>>(
                 theta_steps, phi_steps, pos)) {

            // Loop over volumes
            for (const auto &v : d.volumes()) {
                // Loop over all surfaces in volume
                for (const auto &sf : d.surfaces(v)) {

                    auto sfi = masks.template visit<intersection_update>(
                        sf.mask(), detail::ray(track), sf, transforms);

                    benchmark::DoNotOptimize(hits);
                    benchmark::DoNotOptimize(missed);
                    if (sfi.status == intersection::status::e_inside) {
                        ++hits;
                    } else {
                        ++missed;
                    }
                }
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << "[detray] hits / missed / total = " << hits << " / " << missed
              << " / " << hits + missed << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_ALL)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
