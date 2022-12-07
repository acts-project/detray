/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/geometry/detector_volume.hpp"

/// @note __plugin has to be defined with a preprocessor command

// TODO: Move these into the test defs
namespace {

// geo object ids for testing
enum geo_objects : unsigned int {
    e_sensitive = 0,
    e_portal = 1,
    e_size = 2,
    e_all = e_size,
};

// surface finder ids for testing
enum sf_finder_ids : unsigned int {
    e_default = 0,
};

}  // namespace

// This tests the detector volume class and its many links
TEST(ALGEBRA_PLUGIN, detector_volume) {
    using namespace detray;
    using namespace __plugin;

    using object_link_t = dmulti_index<dindex_range, 2>;
    using sf_finder_link_t = dtyped_index<sf_finder_ids, dindex>;
    using volume_t =
        detector_volume<geo_objects, object_link_t, sf_finder_link_t>;

    // Check construction, setters and getters
    darray<scalar, 6> bounds = {0., 10., -5., 5., -M_PI, M_PI};
    volume_t v1(volume_id::e_cylinder, bounds);
    v1.set_index(12345);
    v1.set_sf_finder({sf_finder_ids::e_default, 12});

    ASSERT_TRUE(v1.empty());
    ASSERT_TRUE(v1.id() == volume_id::e_cylinder);
    ASSERT_TRUE(v1.index() == 12345);
    ASSERT_TRUE(v1.bounds() == bounds);
    ASSERT_TRUE(v1.sf_finder_type() == sf_finder_ids::e_default);
    ASSERT_TRUE(v1.sf_finder_index() == 12);

    // Check surface and portal ranges
    dindex_range surface_range{2, 8};
    dindex_range surface_range_update{8, 10};
    dindex_range full_surface_range{2, 10};
    dindex_range portal_range{10, 24};
    dindex_range full_range{2, 24};

    typename volume_t::sf_finder_link_type sf_finder_link{
        sf_finder_ids::e_default, 12};
    v1.update_obj_link<geo_objects::e_sensitive>(surface_range);
    ASSERT_EQ(v1.obj_link<geo_objects::e_sensitive>(), surface_range);
    v1.update_obj_link<geo_objects::e_sensitive>(surface_range_update);
    ASSERT_EQ(v1.obj_link<geo_objects::e_sensitive>(), full_surface_range);
    v1.update_obj_link<geo_objects::e_portal>(portal_range);
    ASSERT_EQ(v1.obj_link<geo_objects::e_portal>(), portal_range);
    ASSERT_EQ(v1.full_range(), full_range);

    ASSERT_FALSE(v1.empty());
    ASSERT_EQ(v1.template n_objects<geo_objects::e_all>(), 22);
    ASSERT_EQ(v1.template n_objects<geo_objects::e_sensitive>(), 8);
    ASSERT_EQ(v1.template n_objects<geo_objects::e_portal>(), 14);
    ASSERT_EQ(v1.sf_finder_link(), sf_finder_link);

    // Check copy constructor
    const auto v2 = volume_t(v1);
    ASSERT_FALSE(v1.empty());
    ASSERT_EQ(v2.index(), 12345);
    ASSERT_EQ(v2.id(), volume_id::e_cylinder);
    ASSERT_EQ(v2.bounds(), bounds);
    ASSERT_EQ(v2.sf_finder_link(), sf_finder_link);
    ASSERT_EQ(v2.template obj_link<geo_objects::e_sensitive>(),
              full_surface_range);
    ASSERT_EQ(v2.template obj_link<geo_objects::e_portal>(), portal_range);
    ASSERT_EQ(v2.full_range(), full_range);
}
