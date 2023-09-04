/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray core include(s).
#include "detray/masks/masks.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/intersection/intersection.hpp"

// Detray test include(s).
#include "detray/test/types.hpp"

// Google benchmark include(s).
#include <benchmark/benchmark.h>

// System include(s).
#include <iostream>

// Use the detray:: namespace implicitly.
using namespace detray;
using point3 = __plugin::point3<scalar>;

static constexpr unsigned int steps_x3{1000u};
static constexpr unsigned int steps_y3{1000u};
static constexpr unsigned int steps_z3{1000u};

static const test::transform3 trf{};

// This runs a benchmark on a rectangle2D mask
void BM_RECTANGLE_2D_MASK(benchmark::State &state) {

    using mask_type = mask<rectangle2D<>>;
    constexpr mask_type r(0u, 3.f, 4.f);

    constexpr scalar world{10.f};

    constexpr scalar sx{world / steps_x3};
    constexpr scalar sy{world / steps_y3};
    constexpr scalar sz{world / steps_z3};

    unsigned long inside = 0u;
    unsigned long outside = 0u;

    for (auto _ : state) {
        for (unsigned int ix = 0u; ix < steps_x3; ++ix) {
            scalar x{-0.5f * world + static_cast<scalar>(ix) * sx};
            for (unsigned int iy = 0u; iy < steps_y3; ++iy) {
                scalar y{-0.5f * world + static_cast<scalar>(iy) * sy};
                for (unsigned int iz = 0u; iz < steps_z3; ++iz) {
                    scalar z{-0.5f * world + static_cast<scalar>(iz) * sz};

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    const point3 loc_p{r.to_local_frame(trf, {x, y, z})};
                    r.is_inside(loc_p) ? ++inside : ++outside;
                }
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    constexpr scalar area{4.f * r[0] * r[1]};
    constexpr scalar rest{world * world - area};
    std::cout << "Rectangle : Inside/outside ... " << inside << " / " << outside
              << " = "
              << static_cast<scalar>(inside) / static_cast<scalar>(outside)
              << " (theoretical = " << area / rest << ") " << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_RECTANGLE_2D_MASK)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

// This runs a benchmark on a trapezoid2D mask
void BM_TRAPEZOID_2D_MASK(benchmark::State &state) {

    using mask_type = mask<trapezoid2D<>>;
    constexpr mask_type t{0u, 2.f, 3.f, 4.f};

    constexpr scalar world{10.f};

    constexpr scalar sx{world / steps_x3};
    constexpr scalar sy{world / steps_y3};
    constexpr scalar sz{world / steps_z3};

    unsigned long inside = 0u;
    unsigned long outside = 0u;

    for (auto _ : state) {
        for (unsigned int ix = 0u; ix < steps_x3; ++ix) {
            scalar x{-0.5f * world + static_cast<scalar>(ix) * sx};
            for (unsigned int iy = 0u; iy < steps_y3; ++iy) {
                scalar y{-0.5f * world + static_cast<scalar>(iy) * sy};
                for (unsigned int iz = 0u; iz < steps_z3; ++iz) {
                    scalar z{-0.5f * world + static_cast<scalar>(iz) * sz};

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    const point3 loc_p{t.to_local_frame(trf, {x, y, z})};
                    t.is_inside(loc_p) ? ++inside : ++outside;
                }
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    constexpr scalar area{2.f * (t[0] + t[1]) * t[2]};
    constexpr scalar rest{world * world - area};
    std::cout << "Trapezoid : Inside/outside ..." << inside << " / " << outside
              << " = "
              << static_cast<scalar>(inside) / static_cast<scalar>(outside)
              << " (theoretical = " << area / rest << ") " << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_TRAPEZOID_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

// This runs a benchmark on a ring2D mask (as disc)
void BM_DISC_2D_MASK(benchmark::State &state) {

    using mask_type = mask<ring2D<>>;
    constexpr mask_type r{0u, 0.f, 5.f};

    constexpr scalar world{10.f};

    constexpr scalar sx{world / steps_x3};
    constexpr scalar sy{world / steps_y3};
    constexpr scalar sz{world / steps_z3};

    unsigned long inside = 0u;
    unsigned long outside = 0u;

    for (auto _ : state) {
        for (unsigned int ix = 0u; ix < steps_x3; ++ix) {
            scalar x{-0.5f * world + static_cast<scalar>(ix) * sx};
            for (unsigned int iy = 0u; iy < steps_y3; ++iy) {
                scalar y{-0.5f * world + static_cast<scalar>(iy) * sy};
                for (unsigned int iz = 0u; iz < steps_z3; ++iz) {
                    scalar z{-0.5f * world + static_cast<scalar>(iz) * sz};

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    const point3 loc_p{r.to_local_frame(trf, {x, y, z})};
                    r.is_inside(loc_p) ? ++inside : ++outside;
                }
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    constexpr scalar area{r[1] * r[1] * constant<scalar>::pi};
    constexpr scalar rest{world * world - area};
    std::cout << "Disc : Inside/outside ..." << inside << " / " << outside
              << " = "
              << static_cast<scalar>(inside) / static_cast<scalar>(outside)
              << " (theoretical = " << area / rest << ") " << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_DISC_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

// This runs a benchmark on a ring2D mask
void BM_RING_2D_MASK(benchmark::State &state) {

    using mask_type = mask<ring2D<>>;
    constexpr mask_type r{0u, 2.f, 5.f};

    constexpr scalar world{10.f};

    constexpr scalar sx{world / steps_x3};
    constexpr scalar sy{world / steps_y3};
    constexpr scalar sz{world / steps_z3};

    unsigned long inside = 0u;
    unsigned long outside = 0u;

    for (auto _ : state) {
        for (unsigned int ix = 0u; ix < steps_x3; ++ix) {
            scalar x{-0.5f * world + static_cast<scalar>(ix) * sx};
            for (unsigned int iy = 0u; iy < steps_y3; ++iy) {
                scalar y{-0.5f * world + static_cast<scalar>(iy) * sy};
                for (unsigned int iz = 0u; iz < steps_z3; ++iz) {
                    scalar z{-0.5f * world + static_cast<scalar>(iz) * sz};

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    const point3 loc_p{r.to_local_frame(trf, {x, y, z})};
                    r.is_inside(loc_p) ? ++inside : ++outside;
                }
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    constexpr scalar area{(r[1] * r[1] - r[0] * r[0]) * constant<scalar>::pi};
    constexpr scalar rest{world * world - area};
    std::cout << "Ring : Inside/outside ..." << inside << " / " << outside
              << " = "
              << static_cast<scalar>(inside) / static_cast<scalar>(outside)
              << " (theoretical = " << area / rest << ") " << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_RING_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

// This runs a benchmark on a cylinder2D mask
void BM_CYLINDER_2D_MASK(benchmark::State &state) {

    using mask_type = mask<cylinder2D<>>;
    constexpr mask_type c{0u, 3.f, 5.f, 0.f};

    constexpr scalar world{10.f};

    constexpr scalar sx{world / steps_x3};
    constexpr scalar sy{world / steps_y3};
    constexpr scalar sz{world / steps_z3};

    unsigned long inside = 0u;
    unsigned long outside = 0u;

    for (auto _ : state) {
        for (unsigned int ix = 0u; ix < steps_x3; ++ix) {
            scalar x{-0.5f * world + static_cast<scalar>(ix) * sx};
            for (unsigned int iy = 0u; iy < steps_y3; ++iy) {
                scalar y{-0.5f * world + static_cast<scalar>(iy) * sy};
                for (unsigned int iz = 0u; iz < steps_z3; ++iz) {
                    scalar z{-0.5f * world + static_cast<scalar>(iz) * sz};

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    const point3 loc_p{c.to_local_frame(trf, {x, y, z})};
                    c.is_inside(loc_p) ? ++inside : ++outside;
                }
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << "Cylinder : Inside/outside ..." << inside << " / " << outside
              << " = "
              << static_cast<scalar>(inside) / static_cast<scalar>(outside)
              << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_CYLINDER_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

// This runs a benchmark on an annulus2D mask
void BM_ANNULUS_2D_MASK(benchmark::State &state) {

    using mask_type = mask<annulus2D<>>;
    constexpr mask_type ann{0u, 2.5f, 5.f, -0.64299f, 4.13173f, 1.f, 0.5f, 0.f};

    constexpr scalar world{10.f};

    constexpr scalar sx{world / steps_x3};
    constexpr scalar sy{world / steps_y3};
    constexpr scalar sz{world / steps_z3};

    unsigned long inside = 0u;
    unsigned long outside = 0u;

    for (auto _ : state) {
        for (unsigned int ix = 0u; ix < steps_x3; ++ix) {
            scalar x{-0.5f * world + static_cast<scalar>(ix) * sx};
            for (unsigned int iy = 0u; iy < steps_y3; ++iy) {
                scalar y{-0.5f * world + static_cast<scalar>(iy) * sy};
                for (unsigned int iz = 0u; iz < steps_z3; ++iz) {
                    scalar z{-0.5f * world + static_cast<scalar>(iz) * sz};

                    benchmark::DoNotOptimize(inside);
                    benchmark::DoNotOptimize(outside);
                    const point3 loc_p{ann.to_local_frame(trf, {x, y, z})};
                    ann.is_inside(loc_p) ? ++inside : ++outside;
                }
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << "Annulus : Inside/outside ..." << inside << " / " << outside
              << " = "
              << static_cast<scalar>(inside) / static_cast<scalar>(outside)
              << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_ANNULUS_2D_MASK)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
