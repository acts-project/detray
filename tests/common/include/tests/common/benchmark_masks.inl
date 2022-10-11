/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <benchmark/benchmark.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include "detray/intersection/intersection.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

#ifdef DETRAY_BENCHMARKS_REP
unsigned int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
unsigned int gbench_repetitions = 0;
#endif

unsigned int steps_x3 = 1000;
unsigned int steps_y3 = 1000;
unsigned int steps_z3 = 1000;

unsigned int steps_x2 = sqrt(steps_x3 * steps_y3 * steps_z3);
unsigned int steps_y2 = steps_x2;

bool screen_output = false;

using transform_t = __plugin::transform3<detray::scalar>;
using point3_t = __plugin::point3<detray::scalar>;
const transform_t trf{};

namespace {

// This runs a benchmark on a rectangle2D mask
static void BM_RECTANGLE_2D_MASK(benchmark::State &state) {
    using point_t = typename mask<rectangle2D<>>::loc_point_t;

    mask<rectangle2D<>, dindex, transform_t> r{0UL, 3.f, 4.f};

    scalar world = 10.;
    scalar area = 4 * r[0] * r[1];
    scalar rest = world * world - area;

    scalar sx = world / steps_x3;
    scalar sy = world / steps_y3;
    scalar sz = world / steps_z3;

    unsigned long inside = 0;
    unsigned long outside = 0;

    for (auto _ : state) {
        for (unsigned int ix = 0; ix < steps_x3; ++ix) {
            scalar x = -0.5 * world + ix * sx;
            for (unsigned int iy = 0; iy < steps_y3; ++iy) {
                scalar y = -0.5 * world + iy * sy;
                for (unsigned int iz = 0; iz < steps_z3; ++iz) {
                    scalar z = -0.5 * world + iz * sz;

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    point_t loc_p{r.to_local_frame(trf, {x, y, z})};
                    if (r.is_inside(loc_p) == intersection::status::e_inside) {
                        ++inside;
                    } else {
                        ++outside;
                    }
                }
            }
        }
    }

    if (screen_output) {
        std::cout << "Rectangle : Inside/outside ... " << inside << " / "
                  << outside << " = "
                  << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                  << " (theoretical = " << area / rest << ") " << std::endl;
    }
}

// This runs a benchmark on a trapezoid2D mask
static void BM_TRAPEZOID_2D_MASK(benchmark::State &state) {
    using point_t = typename mask<trapezoid2D<>>::loc_point_t;

    mask<trapezoid2D<>> t{0UL, 2.f, 3.f, 4.f};

    scalar world = 10.;
    scalar area = 2 * (t[0] + t[1]) * t[2];
    scalar rest = world * world - area;

    scalar sx = world / steps_x3;
    scalar sy = world / steps_y3;
    scalar sz = world / steps_z3;

    unsigned long inside = 0;
    unsigned long outside = 0;

    for (auto _ : state) {
        for (unsigned int ix = 0; ix < steps_x3; ++ix) {
            scalar x = -0.5 * world + ix * sx;
            for (unsigned int iy = 0; iy < steps_y3; ++iy) {
                scalar y = -0.5 * world + iy * sy;
                for (unsigned int iz = 0; iz < steps_z3; ++iz) {
                    scalar z = -0.5 * world + iz * sz;

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    point_t loc_p{t.to_local_frame(trf, {x, y, z})};
                    if (t.is_inside(loc_p) == intersection::status::e_inside) {
                        ++inside;
                    } else {
                        ++outside;
                    }
                }
            }
        }
    }

    if (screen_output) {
        std::cout << "Trapezoid : Inside/outside ..." << inside << " / "
                  << outside << " = "
                  << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                  << " (theoretical = " << area / rest << ") " << std::endl;
    }
}

// This runs a benchmark on a ring2D mask (as disc)
static void BM_RING_2D_MASK(benchmark::State &state) {
    using point_t = typename mask<ring2D<>>::loc_point_t;

    mask<ring2D<>> r{0UL, 0.f, 5.f};

    scalar world = 10.;
    scalar area = r[1] * r[1] * M_PI;
    scalar rest = world * world - area;

    scalar sx = world / steps_x3;
    scalar sy = world / steps_y3;
    scalar sz = world / steps_z3;

    unsigned long inside = 0;
    unsigned long outside = 0;

    for (auto _ : state) {
        for (unsigned int ix = 0; ix < steps_x3; ++ix) {
            scalar x = -0.5 * world + ix * sx;
            for (unsigned int iy = 0; iy < steps_y3; ++iy) {
                scalar y = -0.5 * world + iy * sy;
                for (unsigned int iz = 0; iz < steps_z3; ++iz) {
                    scalar z = -0.5 * world + iz * sz;

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    point_t loc_p{r.to_local_frame(trf, {x, y, z})};
                    if (r.is_inside(loc_p) == intersection::status::e_inside) {
                        ++inside;
                    } else {
                        ++outside;
                    }
                }
            }
        }
    }
    if (screen_output) {
        std::cout << "Disc : Inside/outside ..." << inside << " / " << outside
                  << " = "
                  << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                  << " (theoretical = " << area / rest << ") " << std::endl;
    }
}

// This runs a benchmark on a ring2D mask
static void BM_DISC_2D_MASK(benchmark::State &state) {
    using point_t = typename mask<ring2D<>>::loc_point_t;

    mask<ring2D<>> r{0UL, 2.f, 5.f};

    scalar world = 10.;
    scalar area = (r[1] * r[1] - r[0] * r[0]) * M_PI;
    scalar rest = world * world - area;

    scalar sx = world / steps_x3;
    scalar sy = world / steps_y3;
    scalar sz = world / steps_z3;

    unsigned long inside = 0;
    unsigned long outside = 0;

    for (auto _ : state) {
        for (unsigned int ix = 0; ix < steps_x3; ++ix) {
            scalar x = -0.5 * world + ix * sx;
            for (unsigned int iy = 0; iy < steps_y3; ++iy) {
                scalar y = -0.5 * world + iy * sy;
                for (unsigned int iz = 0; iz < steps_z3; ++iz) {
                    scalar z = -0.5 * world + iz * sz;

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    point_t loc_p{r.to_local_frame(trf, {x, y, z})};
                    if (r.is_inside(loc_p) == intersection::status::e_inside) {
                        ++inside;
                    } else {
                        ++outside;
                    }
                }
            }
        }
    }
    if (screen_output) {
        std::cout << "Ring : Inside/outside ..." << inside << " / " << outside
                  << " = "
                  << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                  << " (theoretical = " << area / rest << ") " << std::endl;
    }
}

// This runs a benchmark on a cylinder2D mask
static void BM_CYLINDER_2D_MASK(benchmark::State &state) {
    using point_t = typename mask<cylinder2D<>>::loc_point_t;

    mask<cylinder2D<>> c{0UL, 3.f, 5.f, 0.f};

    scalar world = 10.;

    scalar sx = world / steps_x3;
    scalar sy = world / steps_y3;
    scalar sz = world / steps_z3;

    unsigned long inside = 0;
    unsigned long outside = 0;

    for (auto _ : state) {
        for (unsigned int ix = 0; ix < steps_x3; ++ix) {
            scalar x = -0.5 * world + ix * sx;
            for (unsigned int iy = 0; iy < steps_y3; ++iy) {
                scalar y = -0.5 * world + iy * sy;
                for (unsigned int iz = 0; iz < steps_z3; ++iz) {
                    scalar z = -0.5 * world + iz * sz;

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    point_t loc_p{c.to_local_frame(trf, {x, y, z})};
                    if (c.is_inside(loc_p, 0.1) ==
                        intersection::status::e_inside) {
                        ++inside;
                    } else {
                        ++outside;
                    }
                }
            }
        }
    }
    if (screen_output) {
        std::cout << "Cylinder : Inside/outside ..." << inside << " / "
                  << outside << " = "
                  << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                  << std::endl;
    }
}

// This runs a benchmark on an annulus2D mask
static void BM_ANNULUS_2D_MASK(benchmark::State &state) {
    using point_t = typename mask<annulus2D<>>::loc_point_t;

    mask<annulus2D<>> ann{0UL, 2.5f, 5.f, -0.64299f, 4.13173f, 1.f, 0.5f, 0.f};

    scalar world = 10.;

    scalar sx = world / steps_x3;
    scalar sy = world / steps_y3;
    scalar sz = world / steps_z3;

    unsigned long inside = 0;
    unsigned long outside = 0;

    for (auto _ : state) {
        for (unsigned int ix = 0; ix < steps_x3; ++ix) {
            scalar x = -0.5 * world + ix * sx;
            for (unsigned int iy = 0; iy < steps_y3; ++iy) {
                scalar y = -0.5 * world + iy * sy;
                for (unsigned int iz = 0; iz < steps_z3; ++iz) {
                    scalar z = -0.5 * world + iz * sz;

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    point_t loc_p{ann.to_local_frame(trf, {x, y, z})};
                    if (ann.is_inside(loc_p) ==
                        intersection::status::e_inside) {
                        ++inside;
                    } else {
                        ++outside;
                    }
                }
            }
        }
    }
    if (screen_output) {
        std::cout << "Annulus : Inside/outside ..." << inside << " / "
                  << outside << " = "
                  << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                  << std::endl;
    }
}

BENCHMARK(BM_RECTANGLE_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_TRAPEZOID_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_DISC_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_RING_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_CYLINDER_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_ANNULUS_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
#endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

}  // namespace

BENCHMARK_MAIN();
