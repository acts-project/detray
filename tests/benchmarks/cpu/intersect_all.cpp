/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/surface.hpp"
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

using trk_generator_t =
    uniform_track_generator<free_track_parameters<test::transform3>>;

constexpr unsigned int theta_steps{100u};
constexpr unsigned int phi_steps{100u};

// This test runs intersection with all surfaces of the TrackML detector
void BM_INTERSECT_ALL(benchmark::State &state) {

    // Detector configuration
    vecmem::host_memory_resource host_mr;
    toy_det_config toy_cfg{};
    toy_cfg.n_edc_layers(7u);
    auto [d, names] = create_toy_geometry(host_mr, toy_cfg);

    using detector_t = decltype(d);
    detector_t::geometry_context geo_context;

    const auto &transforms = d.transform_store(geo_context);

    std::size_t hits{0u};
    std::size_t missed{0u};
    std::size_t n_surfaces{0u};
    test::point3 origin{0.f, 0.f, 0.f};
    std::vector<intersection2D<typename detector_t::surface_type,
                               typename detector_t::transform3>>
        intersections{};

    // Iterate through uniformly distributed momentum directions
    auto trk_generator = trk_generator_t{};
    trk_generator.config()
        .theta_steps(theta_steps)
        .phi_steps(phi_steps)
        .origin(origin);

    for (auto _ : state) {

        for (const auto track : trk_generator) {

            // Loop over all surfaces in detector
            for (const auto &sf_desc : d.surface_lookup()) {
                const auto sf = surface{d, sf_desc};
                sf.template visit_mask<intersection_initialize>(
                    intersections, detail::ray(track), sf_desc, transforms);

                ++n_surfaces;
            }
            benchmark::DoNotOptimize(hits);
            benchmark::DoNotOptimize(missed);

            hits += intersections.size();
            missed += n_surfaces - intersections.size();

            n_surfaces = 0u;
            intersections.clear();
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
