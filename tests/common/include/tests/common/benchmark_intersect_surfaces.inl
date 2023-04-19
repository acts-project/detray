/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/intersection/concentric_cylinder_intersector.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "tests/common/tools/test_surfaces.hpp"

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// System include(s)
#include <fstream>

using namespace detray;

#ifdef DETRAY_BENCHMARKS_REP
int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
int gbench_repetitions = 0;
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

unsigned int theta_steps = 1000u;
unsigned int phi_steps = 1000u;

dvector<scalar> dists = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f};

/// This benchmark runs intersection with the planar intersector
static void BM_INTERSECT_PLANES(benchmark::State &state) {

    unsigned int sfhit = 0u;
    unsigned int sfmiss = 0u;

    auto planes = planes_along_direction(
        dists, vector::normalize(vector3{1.f, 1.f, 1.f}));
    constexpr mask<rectangle2D<>> rect{0u, 10.f, 20.f};
    point3 ori = {0.f, 0.f, 0.f};

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : uniform_track_generator<detail::ray<transform3>>(
                 theta_steps, phi_steps, ori, 1.f)) {

            for (const auto &plane : planes) {
                auto pi = rect.intersector();
                auto is = pi(ray, plane, rect, plane.transform());

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

/// This benchmark runs intersection with the cylinder intersector
static void BM_INTERSECT_CYLINDERS(benchmark::State &state) {

    using cylinder_mask = mask<cylinder2D<true, cylinder_intersector>>;

    unsigned int sfhit = 0u;
    unsigned int sfmiss = 0u;
    dvector<cylinder_mask> cylinders;

    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask{0u, r, -10.f, 10.f});
    }

    mask_link_t mask_link{mask_ids::e_cylinder2, 0};
    material_link_t material_link{material_ids::e_slab, 0};
    plane_surface plane(transform3(), mask_link, material_link, 0u, false,
                        surface_id::e_sensitive);

    const point3 ori = {0.f, 0.f, 0.f};

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : uniform_track_generator<detail::ray<transform3>>(
                 theta_steps, phi_steps, ori, 1.f)) {

            for (const auto &cylinder : cylinders) {
                auto ci = cylinder.intersector();
                auto inters = ci(ray, plane, cylinder, plane.transform());

                benchmark::DoNotOptimize(sfhit);
                benchmark::DoNotOptimize(sfmiss);
                for (const auto &sfi : inters) {
                    if (sfi.status == intersection::status::e_inside) {
                        ++sfhit;
                    } else {
                        ++sfmiss;
                    }
                }
                benchmark::ClobberMemory();
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

/// This benchmark runs intersection with a specialized portal cylinder
/// intersector
static void BM_INTERSECT_PORTAL_CYLINDERS(benchmark::State &state) {

    using cylinder_mask = mask<cylinder2D<false, cylinder_portal_intersector>>;

    unsigned int sfhit = 0;
    unsigned int sfmiss = 0;
    dvector<cylinder_mask> cylinders;

    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask{0UL, r, -10.f, 10.f});
    }

    mask_link_t mask_link{mask_ids::e_cylinder2, 0};
    material_link_t material_link{material_ids::e_slab, 0};
    plane_surface plane(transform3(), mask_link, material_link, 0, false,
                        surface_id::e_sensitive);

    const point3 ori = {0., 0., 0.};

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : uniform_track_generator<detail::ray<transform3>>(
                 theta_steps, phi_steps, ori, 1.f)) {

            for (const auto &cylinder : cylinders) {
                auto cpi = cylinder.intersector();
                auto is = cpi(ray, plane, cylinder, plane.transform());

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

BENCHMARK(BM_INTERSECT_PORTAL_CYLINDERS)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

/// This benchmark runs intersection with the concentric cylinder intersector
static void BM_INTERSECT_CONCETRIC_CYLINDERS(benchmark::State &state) {
    unsigned int sfhit = 0u;
    unsigned int sfmiss = 0u;

    using cylinder_mask =
        mask<cylinder2D<false, concentric_cylinder_intersector>>;

    dvector<cylinder_mask> cylinders;
    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask(0u, r, -10.f, 10.f));
    }

    mask_link_t mask_link{mask_ids::e_conc_cylinder3, 0u};
    material_link_t material_link{material_ids::e_slab, 0u};
    plane_surface plane(transform3(), mask_link, material_link, 0u, false,
                        surface_id::e_sensitive);

    const point3 ori = {0.f, 0.f, 0.f};

    for (auto _ : state) {

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : uniform_track_generator<detail::ray<transform3>>(
                 theta_steps, phi_steps, ori, 1.f)) {

            for (const auto &cylinder : cylinders) {
                auto cci = cylinder.intersector();
                auto is = cci(ray, plane, cylinder, plane.transform());

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

BENCHMARK(BM_INTERSECT_CONCETRIC_CYLINDERS)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

BENCHMARK_MAIN();
