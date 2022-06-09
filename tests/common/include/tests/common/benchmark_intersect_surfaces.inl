/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <fstream>

#include "detray/core/type_registry.hpp"
#include "detray/intersection/concentric_cylinder_intersector.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/planar_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/test_surfaces.hpp"
#include "tests/common/tools/track_generators.hpp"

using namespace detray;

#ifdef DETRAY_BENCHMARKS_REP
unsigned int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
unsigned int gbench_repetitions = 0;
#endif

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
    e_cylinder3 = 1,
    e_conc_cylinder3 = 2,
};

using mask_defs =
    mask_registry<mask_ids, rectangle2<>, cylinder3<cylinder_intersector>,
                  cylinder3<concentric_cylinder_intersector<>>>;
using plane_surface = surface<mask_defs, transform3>;

unsigned int theta_steps = 1000;
unsigned int phi_steps = 1000;

dvector<scalar> dists = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

// This test runs intersection with all surfaces of the TrackML detector
static void BM_INTERSECT_PLANES(benchmark::State &state) {

    unsigned int sfhit = 0;
    unsigned int sfmiss = 0;

    auto planes =
        planes_along_direction(dists, vector::normalize(vector3{1., 1., 1.}));
    rectangle2<> rect{10., 20., 0u};
    point3 ori = {0., 0., 0.};

    for (auto _ : state) {
        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);

        // Iterate through uniformly distributed momentum directions
        for (const auto test_ray :
             uniform_track_generator<ray>(theta_steps, phi_steps, ori, 1.)) {

            for (auto plane : planes) {
                auto pi = rect.intersector();
                auto is = pi.intersect(plane.transform(), test_ray.pos(),
                                       test_ray.dir(), rect);

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

    using cylinder_mask =
        typename mask_defs::template get_type<e_cylinder3>::type;

    unsigned int sfhit = 0;
    unsigned int sfmiss = 0;
    dvector<cylinder_mask> cylinders;

    for (auto r : dists) {
        cylinders.push_back(cylinder_mask{r, -10., 10., 0u});
    }

    typename mask_defs::link_type mask_link{e_cylinder3, 0};
    plane_surface plain(transform3(), mask_link, 0, false, false);

    point3 ori = {0., 0., 0.};

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

                vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                            cos_theta};

                for (auto cylinder : cylinders) {
                    auto ci = cylinder.intersector();
                    auto is =
                        ci.intersect(plain.transform(), ori, dir, cylinder);

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
        typename mask_defs::template get_type<e_conc_cylinder3>::type;
    ;
    dvector<cylinder_mask> cylinders;
    for (auto r : dists) {
        cylinders.push_back(cylinder_mask{r, -10., 10., 0u});
    }

    typename mask_defs::link_type mask_link{e_conc_cylinder3, 0};
    plane_surface plain(transform3(), mask_link, 0, false, false);

    point3 ori = {0., 0., 0.};

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

                vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                            cos_theta};

                for (auto cylinder : cylinders) {
                    auto cci = cylinder.intersector();
                    auto is =
                        cci.intersect(plain.transform(), ori, dir, cylinder);

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
