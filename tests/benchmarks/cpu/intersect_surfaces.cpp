/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray core include(s).
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/masks/masks.hpp"
#include "detray/navigation/detail/trajectories.hpp"
#include "detray/navigation/intersection/ray_concentric_cylinder_intersector.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"

// Detray utility include(s).
#include "detray/simulation/event_generator/track_generators.hpp"

// Detray test include(s).
#include "detray/test/utils/planes_along_direction.hpp"

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// Use the detray:: namespace implicitly.
using namespace detray;

using ray_generator_t = uniform_track_generator<detail::ray<test::transform3>>;

static const unsigned int theta_steps = 100u;
static const unsigned int phi_steps = 100u;

static const dvector<scalar> dists = {1.f, 2.f, 3.f, 4.f, 5.f,
                                      6.f, 7.f, 8.f, 9.f, 10.f};

/// This benchmark runs intersection with the planar intersector
void BM_INTERSECT_PLANES(benchmark::State &state) {

    unsigned int sfhit = 0u;
    unsigned int sfmiss = 0u;

    auto planes = test::planes_along_direction(
        dists, vector::normalize(test::vector3{1.f, 1.f, 1.f}));
    constexpr mask<rectangle2D> rect{0u, 10.f, 20.f};

    // Iterate through uniformly distributed momentum directions
    auto ray_generator = ray_generator_t{};
    ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : ray_generator) {

            for (const auto &plane : planes) {
                auto pi = ray_intersector<test::transform3, rectangle2D>{};
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
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

namespace {

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
    e_cylinder2 = 1,
    e_conc_cylinder3 = 2,
};

enum material_ids : unsigned int {
    e_slab = 0,
};

}  // namespace

using mask_link_t = dtyped_index<mask_ids, dindex>;
using material_link_t = dtyped_index<material_ids, dindex>;

using plane_surface =
    surface_descriptor<mask_link_t, material_link_t, test::transform3>;
using intersection_t = intersection2D<plane_surface, test::transform3>;

/// This benchmark runs intersection with the cylinder intersector
void BM_INTERSECT_CYLINDERS(benchmark::State &state) {

    using cylinder_mask = mask<cylinder2D>;

    unsigned int sfhit = 0u;
    unsigned int sfmiss = 0u;
    dvector<cylinder_mask> cylinders;

    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask{0u, r, -10.f, 10.f});
    }

    mask_link_t mask_link{mask_ids::e_cylinder2, 0};
    material_link_t material_link{material_ids::e_slab, 0};
    plane_surface plane(test::transform3(), mask_link, material_link, 0u,
                        surface_id::e_sensitive);

    // Iterate through uniformly distributed momentum directions
    auto ray_generator = ray_generator_t{};
    ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : ray_generator) {

            for (const auto &cylinder : cylinders) {
                auto ci = ray_intersector<test::transform3, cylinder2D>{};
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
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with a specialized portal cylinder
/// intersector
void BM_INTERSECT_PORTAL_CYLINDERS(benchmark::State &state) {

    using cylinder_mask = mask<concentric_cylinder2D>;

    unsigned int sfhit = 0u;
    unsigned int sfmiss = 0u;
    dvector<cylinder_mask> cylinders;

    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask{0u, r, -10.f, 10.f});
    }

    mask_link_t mask_link{mask_ids::e_cylinder2, 0u};
    material_link_t material_link{material_ids::e_slab, 0u};
    plane_surface plane(test::transform3(), mask_link, material_link, 0u,
                        surface_id::e_sensitive);

    // Iterate through uniformly distributed momentum directions
    auto ray_generator = ray_generator_t{};
    ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : ray_generator) {

            for (const auto &cylinder : cylinders) {
                auto cpi =
                    ray_intersector<test::transform3, concentric_cylinder2D>{};
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
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the concentric cylinder intersector
void BM_INTERSECT_CONCETRIC_CYLINDERS(benchmark::State &state) {
    unsigned int sfhit = 0u;
    unsigned int sfmiss = 0u;

    using cylinder_mask = mask<concentric_cylinder2D>;

    dvector<cylinder_mask> cylinders;
    for (scalar r : dists) {
        cylinders.push_back(cylinder_mask(0u, r, -10.f, 10.f));
    }

    mask_link_t mask_link{mask_ids::e_conc_cylinder3, 0u};
    material_link_t material_link{material_ids::e_slab, 0u};
    plane_surface plane(test::transform3(), mask_link, material_link, 0u,
                        surface_id::e_sensitive);

    // Iterate through uniformly distributed momentum directions
    auto ray_generator = ray_generator_t{};
    ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

    for (auto _ : state) {

        // Iterate through uniformly distributed momentum directions
        for (const auto ray : ray_generator) {

            for (const auto &cylinder : cylinders) {
                auto cci =
                    ray_concentric_cylinder_intersector<test::transform3>{};
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
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
