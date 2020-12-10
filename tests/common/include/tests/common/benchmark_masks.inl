/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "masks/cylinder3.hpp"
#include "masks/rectangle2.hpp"
#include "masks/trapezoid2.hpp"
#include "masks/ring2.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

#include <benchmark/benchmark.h>

using namespace detray;
using namespace __plugin;

unsigned int steps_x3 = 1000;
unsigned int steps_y3 = 1000;
unsigned int steps_z3 = 1000;

unsigned int steps_x2 = sqrt(steps_x3 * steps_y3 * steps_z3);
unsigned int steps_y2 = steps_x2;

bool screen_output = false;

namespace
{
    // This test runs a rectangle2 maks test operation
    static void BM_RECTANGLE2_MASK(benchmark::State &state)
    {

        using point2 = __plugin::cartesian2::point2;

        rectangle2<scalar> r = {3, 4};

        scalar world = 10.;
        scalar area = 4 * r[0] * r[1];
        scalar rest = world * world - area;

        scalar sx = world / steps_x2;
        scalar sy = world / steps_y2;
        unsigned long inside = 0;
        unsigned long outside = 0;

        for (auto _ : state)
        {
            for (unsigned int ix = 0; ix < steps_x2; ++ix)
            {
                scalar x = -0.5 * world + ix * sx;
                for (unsigned int iy = 0; iy < steps_y2; ++iy)
                {
                    scalar y = -0.5 * world + iy * sy;
                    if (r(point2{x, y}) == e_inside)
                    {
                        ++inside;
                    }
                    else
                    {
                        ++outside;
                    }
                }
            }
        }

        if (screen_output)
        {
            std::cout << "Rectangle : Inside/outside ... " << inside << " / " << outside << " = "
                      << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                      << " (theoretical = " << area / rest << ") " << std::endl;
        }
    }

    // This test runs intersection a trapezoid2 mask operation
    static void BM_TRAPEZOID2_MASK(benchmark::State &state)
    {

        using point2 = __plugin::cartesian2::point2;

        trapezoid2<scalar> t = {2, 3, 4};

        scalar world = 10.;
        scalar area = 2 * (t[0] + t[1]) * t[2];
        scalar rest = world * world - area;

        scalar sx = world / steps_x2;
        scalar sy = world / steps_y2;
        unsigned long inside = 0;
        unsigned long outside = 0;

        for (auto _ : state)
        {
            for (unsigned int ix = 0; ix < steps_x2; ++ix)
            {
                scalar x = -0.5 * world + ix * sx;
                for (unsigned int iy = 0; iy < steps_y2; ++iy)
                {
                    scalar y = -0.5 * world + iy * sy;
                    if (t(point2{x, y}) == e_inside)
                    {
                        ++inside;
                    }
                    else
                    {
                        ++outside;
                    }
                }
            }
        }

        if (screen_output)
        {
            std::cout << "Trapezoid : Inside/outside ..." << inside << " / " << outside << " = "
                      << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                      << " (theoretical = " << area / rest << ") " << std::endl;
        }
    }

    // This test runs a ring2 mask operation
    static void BM_RING2_MASK(benchmark::State &state)
    {

        using point2 = __plugin::cartesian2::point2;

        ring2<scalar> r = {0., 5.};

        scalar world = 10.;
        scalar area = r[1] * r[1] * M_PI;
        scalar rest = world * world - area;

        scalar sx = world / steps_x2;
        scalar sy = world / steps_y2;
        unsigned long inside = 0;
        unsigned long outside = 0;

        for (auto _ : state)
        {
            for (unsigned int ix = 0; ix < steps_x2; ++ix)
            {
                scalar x = -0.5 * world + ix * sx;
                for (unsigned int iy = 0; iy < steps_y2; ++iy)
                {
                    scalar y = -0.5 * world + iy * sy;
                    if (r(point2{x, y}) == e_inside)
                    {
                        ++inside;
                    }
                    else
                    {
                        ++outside;
                    }
                }
            }
        }
        if (screen_output)
        {
            std::cout << "Disc : Inside/outside ..." << inside << " / " << outside << " = "
                      << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                      << " (theoretical = " << area / rest << ") " << std::endl;
        }
    }

    // This test runs mask oeration on a disc2
    static void BM_DISC2_MASK(benchmark::State &state)
    {

        using point2 = __plugin::cartesian2::point2;

        ring2<scalar> r = {2., 5.};

        scalar world = 10.;
        scalar area = (r[1] * r[1] - r[0] * r[0]) * M_PI;
        scalar rest = world * world - area;

        scalar sx = world / steps_x2;
        scalar sy = world / steps_y2;
        unsigned long inside = 0;
        unsigned long outside = 0;

        for (auto _ : state)
        {
            for (unsigned int ix = 0; ix < steps_x2; ++ix)
            {
                scalar x = -0.5 * world + ix * sx;
                for (unsigned int iy = 0; iy < steps_y2; ++iy)
                {
                    scalar y = -0.5 * world + iy * sy;
                    if (r(point2{x, y}) == e_inside)
                    {
                        ++inside;
                    }
                    else
                    {
                        ++outside;
                    }
                }
            }
        }
        if (screen_output)
        {
            std::cout << "Ring : Inside/outside ..." << inside << " / " << outside << " = "
                      << static_cast<scalar>(inside) / static_cast<scalar>(outside)
                      << " (theoretical = " << area / rest << ") " << std::endl;
        }
    }

    // This test runs masking operations on a cylinder3 mask
    static void BM_CYLINDER3_MASK(benchmark::State &state)
    {
        using point3 = __plugin::transform3::point3;

        cylinder3<scalar> c = {3., 5.};

        scalar world = 10.;

        scalar sx = world / steps_x3;
        scalar sy = world / steps_y3;
        scalar sz = world / steps_z3;

        unsigned long inside = 0;
        unsigned long outside = 0;

        for (auto _ : state)
        {
            for (unsigned int ix = 0; ix < steps_x3; ++ix)
            {
                scalar x = -0.5 * world + ix * sx;
                for (unsigned int iy = 0; iy < steps_y3; ++iy)
                {
                    scalar y = -0.5 * world + iy * sy;
                    for (unsigned int iz = 0; iz < steps_z3; ++iz)
                    {
                        scalar z = -0.5 * world + iz * sz;

                        if (c(point3{x, y, z}, 0.1, 0.) == e_inside)
                        {
                            ++inside;
                        }
                        else
                        {
                            ++outside;
                        }
                    }
                }
            }
        }
        if (screen_output)
        {
            std::cout << "Ring : Inside/outside ..." << inside << " / " << outside << " = "
                      << static_cast<scalar>(inside) / static_cast<scalar>(outside) << std::endl;
        }
    }

    BENCHMARK(BM_RECTANGLE2_MASK);
    BENCHMARK(BM_TRAPEZOID2_MASK);
    BENCHMARK(BM_DISC2_MASK);
    BENCHMARK(BM_RING2_MASK);
    BENCHMARK(BM_CYLINDER3_MASK);

} // namespace

BENCHMARK_MAIN();
