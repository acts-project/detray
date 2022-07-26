/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/type_registry.hpp"
#include "detray/geometry/volume.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, volume) {
    using namespace detray;
    using namespace __plugin;

    // dummy types for testing
    enum sf_finder_ids : unsigned int {
        e_default = 0,
    };

    using surface_t = dindex;
    using sf_finder_t = dindex;
    using object_defs = object_registry<surface_t>;
    using sf_finder_defs =
        tuple_array_registry<sf_finder_ids, std::index_sequence<1>,
                             sf_finder_t>;
    using volume = volume<object_defs, scalar, sf_finder_defs::link_type>;

    // Check construction, setters and getters
    darray<scalar, 6> bounds = {0., 10., -5., 5., -M_PI, M_PI};
    volume v1 = volume(bounds);
    v1.set_index(12345);
    v1.set_sf_finder({sf_finder_defs::id::e_default, 12});

    ASSERT_TRUE(v1.empty());
    ASSERT_TRUE(v1.index() == 12345);
    ASSERT_TRUE(v1.bounds() == bounds);
    ASSERT_TRUE(v1.sf_finder_type() == sf_finder_defs::id::e_default);
    ASSERT_TRUE(v1.sf_finder_index() == 12);

    // Check surface and portal ranges
    dindex_range surface_range{2, 8};
    dindex_range portal_range{8, 24};
    dindex_range full_range{2, 24};
    sf_finder_defs::link_type sf_finder_link{sf_finder_defs::id::e_default, 12};
    v1.update_range<object_defs::e_surface>(surface_range);
    v1.update_range<object_defs::e_portal>(portal_range);
    ASSERT_TRUE(v1.range<object_defs::e_surface>() == full_range);
    ASSERT_TRUE(v1.range<object_defs::e_portal>() == full_range);
    ASSERT_FALSE(v1.empty());
    ASSERT_EQ(v1.template n_objects<object_defs::e_surface>(), 22);
    ASSERT_EQ(v1.template n_objects<object_defs::e_portal>(), 22);
    ASSERT_TRUE(v1.sf_finder_link() == sf_finder_link);

    // Check copy constructor
    const auto v2 = volume(v1);
    ASSERT_TRUE(v2.index() == 12345);
    ASSERT_TRUE(v2.bounds() == bounds);
    ASSERT_TRUE(v2.sf_finder_link() == sf_finder_link);
    ASSERT_TRUE(v2.template range<object_defs::e_surface>() == full_range);
    ASSERT_TRUE(v2.template range<object_defs::e_portal>() == full_range);
}
