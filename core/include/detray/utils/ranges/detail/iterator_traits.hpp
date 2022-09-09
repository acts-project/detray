/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <cstddef>
#include <iterator>
#include <type_traits>

#include "detray/definitions/qualifiers.hpp"

namespace detray::ranges::detail {

using std::iterator_traits;

/// Specializations of the std::iterator_traits struct for detray types
/*template<>
struct iterator_traits<transform_store_t<container_t, context_t>> {
    typename difference_type = std::size_t;
    typename value_type = __plugin::transform3<detray::scalar>;;
    typename pointer = __plugin::transform3<detray::scalar> *;
    typename reference = __plugin::transform3<detray::scalar> &;
    typename iterator_category = ;
};*/

}  // namespace detray::ranges::detail