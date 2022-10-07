/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/intersection/line_intersector.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;

namespace {

// 50 mm wire with 1 mm radial cell size
const scalar cell_size{1. * unit_constants::mm};
const scalar hz{50. * unit_constants::mm};

}  // anonymous namespace

/// This tests the basic functionality of a line with a radial cross section
TEST(mask, line_radial_cross_sect) {
    using point_t = typename mask<line<>>::loc_point_t;

    const point_t ln_in{0.09, 0.5};
    const point_t ln_edge{1., 50.};
    const point_t ln_out1{1.2, 0};
    const point_t ln_out2{0.09, -51.};

    const mask<line<>> ln{0UL, cell_size, hz};

    /*ASSERT_FLOAT_EQ(ln[0], static_cast<scalar>(1. * unit_constants::mm));
    ASSERT_FLOAT_EQ(ln[1], static_cast<scalar>(50. * unit_constants::mm));

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_out1) == intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out2) == intersection::status::e_outside);*/
}

/// This tests the basic functionality of a line with a square cross section
TEST(mask, line_square_cross_sect) {
    using point_t = typename mask<line<true>>::loc_point_t;

    /*const point_t ln_in{1., 0., 0};
    const point_t ln_edge{1., 1., 0.};
    const point_t ln_out{1.1, 0., 0};

    // 50 mm wire with 1 mm square cell sizes
    const mask<line<true>> ln{0UL, cell_size, hz};

    ASSERT_FLOAT_EQ(ln[0], static_cast<scalar>(1. * unit_constants::mm));
    ASSERT_FLOAT_EQ(ln[1], static_cast<scalar>(50. * unit_constants::mm));

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, 1e-5) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, -1e-5) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out) == intersection::status::e_outside);*/
}