/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <thrust/tuple.h>

#include <array>
#include <tuple>
#include <utility>

#include "detray/definitions/detray_qualifiers.hpp"

namespace detray {

/** define tuple type namespace for host (std) and device (thrust) compiler
 * **/
#if defined(__CUDACC__)
namespace vtuple = thrust;
#else
namespace vtuple = std;
#endif

namespace detail {

/** tuple size accessor
 */
template <class T>
struct tuple_size;

template <template <typename...> class tuple_type, class... value_types>
struct tuple_size<tuple_type<value_types...>>
    : std::integral_constant<std::size_t, sizeof...(value_types)> {};

/** get function accessor
 *
 */
using std::get;
using thrust::get;

}  // namespace detail
}  // namespace detray