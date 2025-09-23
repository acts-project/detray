/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/type_registry.hpp"

#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <iostream>

namespace detray {

template <bool do_debug = false>
struct select_ray_intersector {
    template <typename mask_t>
    using type = ray_intersector<typename mask_t::shape,
                                 typename mask_t::algebra_type, do_debug>;
};

}  // namespace detray

// Test type list implementation
GTEST_TEST(detray_utils, mapped_type_registry) {
    using namespace detray;

    using metadata_t = test::toy_metadata;
    using detector_t = detector<metadata_t>;
    using test_algebra_t = detector_t::algebra_type;
    using mask_id = typename detector_t::masks::id;

    using mapped_registry_t =
        mapped_type_registry<typename detector_t::masks,
                             select_ray_intersector<true>>;

    types::print<types::list<typename mapped_registry_t::types>>();
    std::cout << std::endl;

    //
    // Test the mapping
    //
    constexpr const auto& id_array = mapped_registry_t::ids();

    static_assert(mapped_registry_t::n_types == 3u);
    static_assert(id_array.size() == 4u);

    EXPECT_EQ(id_array[0u], 0u);
    EXPECT_EQ(id_array[1u], 0u);
    EXPECT_EQ(id_array[2u], 1u);
    EXPECT_EQ(id_array[3u], 2u);

    using cyl_intersector_t = typename mapped_registry_t::template get_type<
        mask_id::e_cylinder2>::type;

    using rect_intersector_t = typename mapped_registry_t::template get_type<
        mask_id::e_rectangle2>::type;

    using ring_intersector_t = typename mapped_registry_t::template get_type<
        mask_id::e_portal_ring2>::type;

    using trpz_intersector_t = typename mapped_registry_t::template get_type<
        mask_id::e_trapezoid2>::type;

    types::print<types::list<cyl_intersector_t>>();
    std::cout << std::endl;
    types::print<types::list<rect_intersector_t>>();
    std::cout << std::endl;
    types::print<types::list<ring_intersector_t>>();
    std::cout << std::endl;
    types::print<types::list<trpz_intersector_t>>();

    static_assert(
        std::same_as<cyl_intersector_t, ray_intersector<concentric_cylinder2D,
                                                        test_algebra_t, true>>,
        "Retrieved incorrect type after mapping (cylinder2D)");

    static_assert(
        std::same_as<rect_intersector_t,
                     ray_intersector<rectangle2D, test_algebra_t, true>>,
        "Retrieved incorrect type after mapping (rectangle2D)");

    static_assert(std::same_as<ring_intersector_t,
                               ray_intersector<ring2D, test_algebra_t, true>>,
                  "Retrieved incorrect type after mapping (ring2D)");

    static_assert(
        std::same_as<trpz_intersector_t,
                     ray_intersector<trapezoid2D, test_algebra_t, true>>,
        "Retrieved incorrect type after mapping (trapezoid2D)");

    //
    // Test the registry
    //

    // Get ID
    static_assert(mapped_registry_t::get_id<cyl_intersector_t>() ==
                      static_cast<mask_id>(
                          id_array[static_cast<dindex>(mask_id::e_cylinder2)]),
                  "ID for type cylinder intersector incorrect");
    static_assert(mapped_registry_t::get_id<rect_intersector_t>() ==
                      static_cast<mask_id>(
                          id_array[static_cast<dindex>(mask_id::e_rectangle2)]),
                  "ID for type rectangle intersector incorrect");
    static_assert(
        mapped_registry_t::get_id<ring_intersector_t>() ==
            static_cast<mask_id>(
                id_array[static_cast<dindex>(mask_id::e_portal_ring2)]),
        "ID for type ring intersector incorrect");
    static_assert(mapped_registry_t::get_id<trpz_intersector_t>() ==
                      static_cast<mask_id>(
                          id_array[static_cast<dindex>(mask_id::e_trapezoid2)]),
                  "ID for type trapezoid intersector incorrect");

    // contains
    static_assert(mapped_registry_t::template contains<cyl_intersector_t>(),
                  "'contains' failed for cylinder intersector type");
    static_assert(mapped_registry_t::template contains<rect_intersector_t>(),
                  "'contains' failed for rectangle intersector type");
    static_assert(mapped_registry_t::template contains<ring_intersector_t>(),
                  "'contains' failed for ring intersector type");
    static_assert(mapped_registry_t::template contains<trpz_intersector_t>(),
                  "'contains' failed for trapezoid intersector type");

    static_assert(!mapped_registry_t::template contains<char>(),
                  "'contains' failed for 'char' type");
    static_assert(!mapped_registry_t::template contains<void>(),
                  "'contains' failed for 'void' type");

    // Is valid
    static_assert(mapped_registry_t::is_valid(mask_id::e_cylinder2),
                  "ID for cylinder intersector invalid");
    static_assert(mapped_registry_t::is_valid(mask_id::e_rectangle2),
                  "ID for rectangle intersector invalid");
    static_assert(mapped_registry_t::is_valid(mask_id::e_portal_ring2),
                  "ID for ring intersector invalid");
    static_assert(mapped_registry_t::is_valid(mask_id::e_trapezoid2),
                  "ID for trapezoid intersector invalid");
    assert(!mapped_registry_t::is_valid(5u) && "ID '5' not invalid");

    // Get type
    static_assert(std::same_as<cyl_intersector_t,
                               typename mapped_registry_t::template get_type<
                                   mask_id::e_cylinder2>::type>,
                  "Got incorrect type for 'e_cylinder2'");
    static_assert(std::same_as<rect_intersector_t,
                               typename mapped_registry_t::template get_type<
                                   mask_id::e_rectangle2>::type>,
                  "Got incorrect type for 'e_rectangle2'");
    static_assert(std::same_as<ring_intersector_t,
                               typename mapped_registry_t::template get_type<
                                   mask_id::e_portal_ring2>::type>,
                  "Got incorrect type for 'e_portal_ring2'");
    static_assert(std::same_as<trpz_intersector_t,
                               typename mapped_registry_t::template get_type<
                                   mask_id::e_trapezoid2>::type>,
                  "Got incorrect type for 'e_trapezoid2'");
}
