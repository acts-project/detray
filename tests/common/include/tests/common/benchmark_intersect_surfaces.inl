/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/intersection/concentric_cylinder_intersector.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/simulation/track_generators.hpp"
#include "tests/common/tools/test_surfaces.hpp"

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// System include(s)
#include <fstream>

using namespace detray;

#ifdef DETRAY_BENCHMARKS_REP
unsigned int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
unsigned int gbench_repetitions = 0;
#endif

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
    e_cylinder2 = 1,
    e_conc_cylinder3 = 2,
};

enum material_ids : unsigned int {
    e_slab = 0,
};

using mask_link_t = dtyped_index<mask_ids, dindex>;
using material_link_t = dtyped_index<material_ids, dindex>;

using plane_surface = surface<mask_link_t, material_link_t, transform3>;

unsigned int theta_steps = 1000;
unsigned int phi_steps = 1000;

dvector<scalar> dists = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

// This test runs intersection with all surfaces of the TrackML detector
static void BM_INTERSECT_PLANES(benchmark::State &state) {

    unsigned int sfhit = 0;
    unsigned int sfmiss = 0;

    auto planes =
        planes_along_direction(dists, vector::normalize(vector3{1., 1., 1.}));
    mask<rectangle2D<>> rect{0UL, 10.f, 20.f};
    point3 ori = {0., 0., 0.};

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : uniform_track_generator<detail::ray<transform3>>(
                 theta_steps, phi_steps, ori, 1.)) {

            for (const auto &plane : planes) {
                auto pi = rect.intersector();
                auto is = pi(ray, rect, plane.transform())[0];

                benchmark::DoNotOptimize(sfhit);
                benchmark::DoNotOptimize(sfmiss);
                if (is.status == intersection::status::e_inside) {
                    ++sfhit;
                } else {
                    ++sfmiss;
                }
                benchmark::ClobberMemory();
            }
        }
    }
}

BENCHMARK(BM_INTERSECT_PLANES)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

// This test runs intersection with all surfaces of the TrackML detector
static void BM_INTERSECT_CYLINDERS(benchmark::State &state) {

    using cylinder_mask = mask<cylinder2D<>>;

    unsigned int sfhit = 0;
    unsigned int sfmiss = 0;
    dvector<cylinder_mask> cylinders;

    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask{0UL, r, -10.f, 10.f});
    }

    mask_link_t mask_link{mask_ids::e_cylinder2, 0};
    material_link_t material_link{material_ids::e_slab, 0};
    plane_surface plain(transform3(), mask_link, material_link, 0, false,
                        surface_id::e_sensitive);

    const point3 ori = {0., 0., 0.};

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Loops of theta values
        for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
            scalar theta = 0.85 + itheta * 0.2 / theta_steps;
            scalar sin_theta = std::sin(theta);
            scalar cos_theta = std::cos(theta);

            // Loops of phi values
            for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
                scalar phi = 0.7 + iphi * 0.2 / phi_steps;
                scalar sin_phi = std::sin(phi);
                scalar cos_phi = std::cos(phi);

                const vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                                  cos_theta};

                const detail::ray<transform3> ray(ori, 0., dir, 0.);

                for (const auto &cylinder : cylinders) {
                    auto ci = cylinder.intersector();
                    auto is = ci(ray, cylinder, plain.transform())[0];

                    benchmark::DoNotOptimize(sfhit);
                    benchmark::DoNotOptimize(sfmiss);
                    if (is.status == intersection::status::e_inside) {
                        ++sfhit;
                    } else {
                        ++sfmiss;
                    }
                    benchmark::ClobberMemory();
                }
            }
        }
    }
}

BENCHMARK(BM_INTERSECT_CYLINDERS)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

// This test runs intersection with all surfaces of the TrackML detector
static void BM_INTERSECT_CONCETRIC_CYLINDERS(benchmark::State &state) {
    unsigned int sfhit = 0;
    unsigned int sfmiss = 0;

    using cylinder_mask =
        mask<cylinder2D<false, concentric_cylinder_intersector>>;

    dvector<cylinder_mask> cylinders;
    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask(0UL, r, -10.f, 10.f));
    }

    mask_link_t mask_link{mask_ids::e_conc_cylinder3, 0};
    material_link_t material_link{material_ids::e_slab, 0};
    plane_surface plain(transform3(), mask_link, material_link, 0, false,
                        surface_id::e_sensitive);

    const point3 ori = {0., 0., 0.};

    for (auto _ : state) {
        // Loops of theta values
        for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
            scalar theta = 0.85 + itheta * 0.2 / theta_steps;
            scalar sin_theta = std::sin(theta);
            scalar cos_theta = std::cos(theta);

            // Loops of phi values
            for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
                scalar phi = 0.7 + iphi * 0.2 / phi_steps;
                scalar sin_phi = std::sin(phi);
                scalar cos_phi = std::cos(phi);

                const vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                                  cos_theta};

                const detail::ray<transform3> ray(ori, 0., dir, 0.);

                for (const auto &cylinder : cylinders) {
                    auto cci = cylinder.intersector();
                    auto is = cci(ray, cylinder, plain.transform())[0];

                    benchmark::DoNotOptimize(sfhit);
                    benchmark::DoNotOptimize(sfmiss);
                    if (is.status == intersection::status::e_inside) {
                        ++sfhit;
                    } else {
                        ++sfmiss;
                    }
                    benchmark::ClobberMemory();
                }
            }
        }
    }
}

BENCHMARK(BM_INTERSECT_CONCETRIC_CYLINDERS)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

BENCHMARK_MAIN();
