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

/** get function accessor
 *
 */
using std::get;
using thrust::get;

/** tuple size accessor
 */
template <class T>
struct tuple_size;

template <template <typename...> class tuple_type, class... value_types>
struct tuple_size<tuple_type<value_types...>>
    : std::integral_constant<std::size_t, sizeof...(value_types)> {};

/** make tuple accessor
 *  users have to specifiy tuple_type for detail::make_tuple
 */

template <class T>
struct unwrap_refwrapper {
    using type = T;
};

template <class T>
struct unwrap_refwrapper<std::reference_wrapper<T>> {
    using type = T&;
};

template <class T>
using unwrap_decay_t =
    typename unwrap_refwrapper<typename std::decay<T>::type>::type;

template <template <typename...> class tuple_type, class... Types>
constexpr tuple_type<unwrap_decay_t<Types>...> make_tuple(Types&&... args) {
    return tuple_type<unwrap_decay_t<Types>...>(std::forward<Types>(args)...);
}

}  // namespace detail
}  // namespace detray