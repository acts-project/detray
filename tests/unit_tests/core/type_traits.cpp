/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"

// Vecmem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/containers/vector.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <type_traits>

TEST(type_traits, is_same_vector) {

    using namespace detray;

    static_assert(is_same_vector_v<std::vector, std::vector> == true);
    static_assert(is_same_vector_v<vecmem::vector, vecmem::vector> == true);
    static_assert(
        is_same_vector_v<vecmem::jagged_vector, vecmem::jagged_vector> == true);
    static_assert(
        is_same_vector_v<vecmem::device_vector, vecmem::device_vector> == true);
    static_assert(is_same_vector_v<vecmem::jagged_device_vector,
                                   vecmem::jagged_device_vector> == true);
    static_assert(is_same_vector_v<std::vector, vecmem::vector> == false);
    static_assert(is_same_vector_v<vecmem::vector, vecmem::device_vector> ==
                  false);
}