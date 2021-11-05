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

/** Tuple of vector which can hold either host and device vector
 *
 */
template <template <typename...> class vector_type = vecmem::vector,
          typename... Ts>
struct vec_tuple {

    // vtuple has different types based on the file location
    // std::tuple in cpp/hpp;
    // thrust::tuple in cu
    vtuple::tuple<vector_type<Ts>...> _tuple;

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
    template <template <typename...> class vec_tuple_data>
    DETRAY_DEVICE vec_tuple(vec_tuple_data<Ts...>& data)
        : _tuple(data.device(std::make_index_sequence<sizeof...(Ts)>{})) {}
};

/** Tuple of vecmem::vector_view
 *
 */
template <typename... Ts>
struct vec_tuple_data {
    thrust::tuple<vecmem::data::vector_view<Ts>...> _tuple;

    /** Constructor with vec_tuple - only for host
     *
     */
    template <template <typename...> class vector_type, std::size_t... ints>
    DETRAY_HOST vec_tuple_data(vec_tuple<vector_type, Ts...>& vtuple,
                               std::index_sequence<ints...> /*seq*/)
        : _tuple(thrust::make_tuple(vecmem::data::vector_view<Ts>(
              vecmem::get_data(std::get<ints>(vtuple._tuple)))...)) {}

    /** Create tuple with vecmem::device_vector
     *
     * This function is called by vec_tuple constructor
     */
    template <std::size_t... ints>
    DETRAY_DEVICE thrust::tuple<vecmem::device_vector<Ts>...> device(
        std::index_sequence<ints...> /*seq*/) {
        return thrust::make_tuple(
            vecmem::device_vector<Ts>(thrust::get<ints>(_tuple))...);
    }
};

/** Get vec_tuple_data
 **/
template <template <typename...> class vector_type, typename... Ts>
inline vec_tuple_data<Ts...> get_data(vec_tuple<vector_type, Ts...>& vtuple) {
    return vec_tuple_data<Ts...>(vtuple,
                                 std::make_index_sequence<sizeof...(Ts)>{});
}

// Test function to copy the contenst of vec_tuple_data into vecmem vector
void tuple_copy(vec_tuple_data<int, float, double>& data,
                vecmem::data::vector_view<int>& output1_data,
                vecmem::data::vector_view<float>& output2_data,
                vecmem::data::vector_view<double>& output3_data);
}  // namespace detray
