/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/type_registry.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <cstdint>
#include <type_traits>

// Test IDs
enum class test_ids : std::uint_least8_t {
    e_int = 0u,
    e_double = 1u,
    e_float1 = 2u,
    e_float2 = 3u,
};

// Test type list implementation
GTEST_TEST(detray_utils, type_registry) {
    using namespace detray;

    using type_registry_t = type_registry<test_ids, int, double, float, float>;

    static_assert(type_registry_t::n_types == 4u,
                  "Incorrect number of types in registry");

    // Get ID
    static_assert(type_registry_t::get_id<int>() == static_cast<test_ids>(0u),
                  "ID for type 'int' incorrect");
    static_assert(
        type_registry_t::get_id<double>() == static_cast<test_ids>(1u),
        "ID for type 'double' incorrect");
    static_assert(type_registry_t::get_id<float>() == static_cast<test_ids>(2u),
                  "ID for type 'float' incorrect");

    // contains
    static_assert(type_registry_t::template contains<int>(),
                  "'contains' failed for 'int' type");
    static_assert(type_registry_t::template contains<double>(),
                  "'contains' failed for 'double' type");
    static_assert(type_registry_t::template contains<float>(),
                  "'contains' failed for 'float' type");

    static_assert(!type_registry_t::template contains<char>(),
                  "'contains' failed for 'char' type");
    static_assert(!type_registry_t::template contains<void>(),
                  "'contains' failed for 'void' type");

    // Is valid
    static_assert(
        type_registry_t::is_valid(static_cast<std::size_t>(test_ids::e_int)),
        "ID for 'int' invalid");
    static_assert(
        type_registry_t::is_valid(static_cast<std::size_t>(test_ids::e_double)),
        "ID for 'double' invalid");
    static_assert(
        type_registry_t::is_valid(static_cast<std::size_t>(test_ids::e_float1)),
        "ID for 'float 1' invalid");
    static_assert(
        type_registry_t::is_valid(static_cast<std::size_t>(test_ids::e_float2)),
        "ID for 'float 2' invalid");
    static_assert(!type_registry_t::is_valid(5u), "ID '5' not invalid");

    // Get type
    static_assert(std::same_as<int, typename type_registry_t::template get_type<
                                        test_ids::e_int>::type>,
                  "Got incorrect type for 'e_int'");
    static_assert(
        std::same_as<double, typename type_registry_t::template get_type<
                                 test_ids::e_double>::type>,
        "Got incorrect type for 'e_double'");
    static_assert(
        std::same_as<float, typename type_registry_t::template get_type<
                                test_ids::e_float1>::type>,
        "Got incorrect type for 'e_float1'");
    static_assert(
        std::same_as<float, typename type_registry_t::template get_type<
                                test_ids::e_float2>::type>,
        "Got incorrect type for 'e_float2'");
}
