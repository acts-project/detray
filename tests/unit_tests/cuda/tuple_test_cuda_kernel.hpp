/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#if defined(array)
#include "detray/plugins/algebra/array_definitions.hpp"
#elif defined(eigen)
#include "detray/plugins/algebra/eigen_definitions.hpp"
#elif defined(smatrix)
#include "detray/plugins/algebra/smatrix_definitions.hpp"
#elif defined(vc_array)
#include "detray/plugins/algebra/vc_array_definitions.hpp"
#endif

#include <vecmem/containers/device_vector.hpp>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/detray_qualifiers.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

namespace detail {

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
// or use std::unwrap_ref_decay_t (since C++20)

template <template <typename...> class tuple_type, class... Types>
constexpr  // since C++14
    tuple_type<unwrap_decay_t<Types>...>
    make_tuple(Types&&... args) {
    return tuple_type<unwrap_decay_t<Types>...>(std::forward<Types>(args)...);
}

}  // namespace detail

/** Tuple of vector which can hold either host and device vector
 *
 */
template <template <typename...> class tuple_type = std::tuple,
          template <typename...> class vector_type = vecmem::vector,
          typename... Ts>
struct vec_tuple {

    // vtuple has different types based on the file location
    // std::tuple in cpp/hpp;
    // thrust::tuple in cu
    vtuple::tuple<vector_type<Ts>...> _tuple;

    // tuple_type for _tuple makes an illegal memory access error
    // tuple_type<vector_type<Ts>...> _tuple;

    /** Constructor with vec_tuple_data - only for device
     *
     * Used when creating vec_tuple with vecmem::vector
     */
    DETRAY_HOST vec_tuple(vecmem::memory_resource& resource)
        : _tuple(vector_type<Ts>{&resource}...) {}

    /** Constructor with vec_tuple_data - only for device
     *
     * Used when creating vec_tuple with vecmem::device_vector
     */
    template <template <template <typename...> class, typename...>
              class vec_tuple_data>
    DETRAY_DEVICE vec_tuple(vec_tuple_data<tuple_type, Ts...>& data)
        : _tuple(data.device(std::make_index_sequence<sizeof...(Ts)>{})) {}
};

/** Tuple of vecmem::vector_view
 *
 */
template <template <typename...> class tuple_type, typename... Ts>
struct vec_tuple_data {
    tuple_type<vecmem::data::vector_view<Ts>...> _tuple;

    /** Constructor with vec_tuple - only for host
     *
     */
    template <template <typename...> class vector_type, std::size_t... ints>
    DETRAY_HOST vec_tuple_data(
        vec_tuple<tuple_type, vector_type, Ts...>& vtuple,
        std::index_sequence<ints...> /*seq*/)
        : _tuple(detail::make_tuple<tuple_type>(vecmem::data::vector_view<Ts>(
              vecmem::get_data(detail::get<ints>(vtuple._tuple)))...)) {}

    /** Create tuple with vecmem::device_vector
     *
     * This function is called by vec_tuple constructor
     */
    template <std::size_t... ints>
    DETRAY_DEVICE tuple_type<vecmem::device_vector<Ts>...> device(
        std::index_sequence<ints...> /*seq*/) {
        return detail::make_tuple<tuple_type>(
            vecmem::device_vector<Ts>(detail::get<ints>(_tuple))...);
    }
};

/** Get vec_tuple_data
 **/
template <template <typename...> class tuple_type,
          template <typename...> class vector_type, typename... Ts>
inline vec_tuple_data<tuple_type, Ts...> get_data(
    vec_tuple<tuple_type, vector_type, Ts...>& vtuple) {
    return vec_tuple_data<tuple_type, Ts...>(
        vtuple, std::make_index_sequence<sizeof...(Ts)>{});
}

// Test function to copy the contenst of vec_tuple_data into vecmem vector
template <template <typename...> class tuple_type>
void tuple_copy(vec_tuple_data<tuple_type, int, float, double>& data,
                vecmem::data::vector_view<int>& output1_data,
                vecmem::data::vector_view<float>& output2_data,
                vecmem::data::vector_view<double>& output3_data);
}  // namespace detray
