/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <vecmem/memory/host_memory_resource.hpp>

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include <gtest/gtest.h>

#include <climits>

#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/utils/indexing.hpp"

using namespace detray;

TEST(grids, grid2_replace_populator) {
    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r = grid2<replace_populator, axis::regular, axis::regular,
                         decltype(serializer)>;
    typename grid2r::axis_p0_type xaxis{10, -5., 5., host_mr};
    typename grid2r::axis_p1_type yaxis{10, -5., 5., host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr);

    // Test the initialization
    test::point2<detray::scalar> p = {-4.5, -4.5};
    for (unsigned int ib0 = 0; ib0 < 10; ++ib0) {
        for (unsigned int ib1 = 0; ib1 < 10; ++ib1) {
            p = {static_cast<scalar>(-4.5 + ib0),
                 static_cast<scalar>(-4.5 + ib1)};
            EXPECT_EQ(g2.bin(p), std::numeric_limits<dindex>::max());
        }
    }

    p = {-4.5, -4.5};
    // Fill and read
    g2.populate(p, 3u);
    EXPECT_EQ(g2.bin(p), 3u);

    // Fill and read two times, fill first 0-99, then 100-199
    for (unsigned int il = 0; il < 2; ++il) {
        unsigned int counter = il * 100;
        for (unsigned int ib0 = 0; ib0 < 10; ++ib0) {
            for (unsigned int ib1 = 0; ib1 < 10; ++ib1) {
                p = {static_cast<scalar>(-4.5 + ib0),
                     static_cast<scalar>(-4.5 + ib1)};
                g2.populate(p, counter);
                EXPECT_EQ(g2.bin(p), counter++);
            }
        }
    }

    // A zone test w/o neighbour hood
    p = {-4.5, -4.5};
    auto test = g2.zone(p);
    dvector<dindex> expect = {100u};
    EXPECT_EQ(test, expect);

    // A zone test with neighbour hood
    p = {0.5, 0.5};

    darray<dindex, 2> zone11 = {1u, 1u};
    darray<dindex, 2> zone22 = {2u, 2u};

    test = g2.zone(p, {zone11, zone22}, true);
    expect = {143u, 144u, 145u, 146u, 147u, 153u, 154u, 155u,
              156u, 157u, 163u, 164u, 165u, 166u, 167u};
    EXPECT_EQ(test, expect);

    using grid2cc = grid2<replace_populator, axis::circular, axis::regular,
                          decltype(serializer)>;

    typename grid2cc::axis_p0_type circular{4, -2., 2., host_mr};
    typename grid2cc::axis_p1_type closed{5, 0., 5., host_mr};

    grid2cc g2cc(std::move(circular), std::move(closed), host_mr);
    unsigned int counter = 0;
    for (unsigned icl = 0; icl < 5; ++icl) {
        for (unsigned ici = 0; ici < 4; ++ici) {
            p = {static_cast<scalar>(-1.5 + ici),
                 static_cast<scalar>(0.5 + icl)};
            g2cc.populate(p, counter++);
        }
    }

    // A zone test for circular testing
    p = {1.5, 2.5};
    test = g2cc.zone(p, {zone11, zone11}, true);
    expect = {4u, 6u, 7u, 8u, 10u, 11u, 12u, 14u, 15u};
    EXPECT_EQ(test, expect);
}

TEST(grids, grid2_complete_populator) {
    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r = grid2<complete_populator, axis::regular, axis::regular,
                         decltype(serializer), dvector, djagged_vector, darray,
                         dtuple, dindex, false, 3>;

    typename grid2r::axis_p0_type xaxis{2, -1., 1., host_mr};
    typename grid2r::axis_p1_type yaxis{2, -1., 1., host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr);

    // Test the initialization
    test::point2<detray::scalar> p = {-0.5, -0.5};
    grid2r::populator_type::store_value invalid = {
        dindex_invalid, dindex_invalid, dindex_invalid};
    for (unsigned int ib0 = 0; ib0 < 2; ++ib0) {
        for (unsigned int ib1 = 0; ib1 < 2; ++ib1) {
            p = {static_cast<scalar>(-0.5 + ib0),
                 static_cast<scalar>(-0.5 + ib1)};
            EXPECT_EQ(g2.bin(p), invalid);
        }
    }

    // Fill and read
    p = {-0.5, -0.5};
    g2.populate(p, 4u);

    grid2r::populator_type::store_value expected = {4u, dindex_invalid,
                                                    dindex_invalid};
    auto test = g2.bin(p);
    EXPECT_EQ(test, expected);

    auto zone_test = g2.zone(p);
    dvector<dindex> zone_expected = {4u};
    EXPECT_EQ(zone_test, zone_expected);

    g2.populate(p, 2u);
    expected = {4u, 2u, dindex_invalid};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    g2.populate(p, 7u);
    expected = {4u, 2u, 7u};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    // Bin is completed, new entry is ignored
    g2.populate(p, 16u);
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    darray<dindex, 2> zone00 = {0u, 0u};
    darray<dindex, 2> zone11 = {1u, 1u};

    // Zone test of a complete bin
    zone_test = g2.zone(p, {zone00, zone00});
    zone_expected = {4u, 2u, 7u};
    EXPECT_EQ(zone_test, zone_expected);

    // Fill some other bins
    p = {0.5, -0.5};
    g2.populate(p, 16u);

    p = {0.5, 0.5};
    g2.populate(p, 17u);
    g2.populate(p, 18u);

    zone_test = g2.zone(p, {zone11, zone11});
    zone_expected = {4u, 2u, 7u, 16u, 17u, 18u};
    EXPECT_EQ(zone_test, zone_expected);
}

TEST(grids, grid2_attach_populator) {
    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r = grid2<attach_populator, axis::regular, axis::regular,
                         decltype(serializer)>;
    typename grid2r::axis_p0_type xaxis{2, -1., 1., host_mr};
    typename grid2r::axis_p1_type yaxis{2, -1., 1., host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr);

    // Test the initialization
    test::point2<detray::scalar> p = {-0.5, -0.5};
    grid2r::populator_type::store_value invalid = {};
    for (unsigned int ib0 = 0; ib0 < 2; ++ib0) {
        for (unsigned int ib1 = 0; ib1 < 2; ++ib1) {
            p = {static_cast<scalar>(-0.5 + ib0),
                 static_cast<scalar>(-0.5 + ib1)};
            EXPECT_EQ(g2.bin(p), invalid);
        }
    }

    p = {-0.5, -0.5};
    g2.populate(p, 4u);

    grid2r::populator_type::store_value expected = {4u};
    auto test = g2.bin(p);
    EXPECT_EQ(test, expected);

    auto zone_test = g2.zone(p);
    dvector<dindex> zone_expected = {4u};
    EXPECT_EQ(zone_test, zone_expected);

    p = {-0.5, 0.5};
    g2.populate(p, 9u);

    p = {0.5, -0.5};
    g2.populate(p, 1u);

    p = {0.5, 0.5};
    g2.populate(p, 7u);

    expected = {7u};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    darray<dindex, 2> zone11 = {1u, 1u};

    zone_test = g2.zone(p, {zone11, zone11}, true);
    zone_expected = {1u, 4u, 7u, 9u};
    EXPECT_EQ(zone_test, zone_expected);
}

TEST(grids, grid2_shift) {
    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r = grid2<replace_populator, axis::regular, axis::regular,
                         decltype(serializer)>;

    typename grid2r::axis_p0_type xaxis{10, -5., 5., host_mr};
    typename grid2r::axis_p1_type yaxis{10, -5., 5., host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr, 0);

    // Test the initialization
    test::point2<detray::scalar> p = {-4.5, -4.5};
    EXPECT_EQ(g2.bin(p), 0u);

    g2.shift(8u);
    EXPECT_EQ(g2.bin(p), 8u);
}

TEST(grids, grid2_irregular_replace) {
    vecmem::host_memory_resource host_mr;

    replace_populator<> replacer;
    serializer2 serializer;

    using grid2ir = grid2<replace_populator, axis::irregular, axis::irregular,
                          decltype(serializer)>;

    typename grid2ir::axis_p0_type xaxis{
        {-3, -2., 1, 0.5, 0.7, 0.71, 4., 1000.}, host_mr};
    typename grid2ir::axis_p1_type yaxis{{0.1, 0.8, 0.9, 10., 12., 15.},
                                         host_mr};

    grid2ir g2(std::move(xaxis), std::move(yaxis), host_mr);

    test::point2<detray::scalar> p = {-0.5, 0.5};
    g2.populate(p, 4u);
    EXPECT_EQ(g2.bin(p), 4u);
}
