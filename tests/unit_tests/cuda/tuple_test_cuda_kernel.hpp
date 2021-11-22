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
#include "detray/definitions/qualifiers.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

/** Tuple of vector which can hold either host and device vector
 *
 */
template <template <typename...> class tuple_type = std::tuple,
          template <typename...> class vector_type = vecmem::vector,
          typename... Ts>
struct vec_tuple {

    using data_t = tuple_type<vecmem::data::vector_view<Ts>...>;

    // vtuple has different types based on the file location
    // 1) std::tuple in cpp/hpp; 2) thrust::tuple in cu
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
    template <typename vec_tuple_data_t,
              std::enable_if_t<
                  !std::is_base_of_v<vecmem::memory_resource, vec_tuple_data_t>,
                  bool> = true>
    DETRAY_DEVICE vec_tuple(vec_tuple_data_t& data)
        : _tuple(device(data, std::make_index_sequence<sizeof...(Ts)>{})) {}

    template <std::size_t... ints>
    DETRAY_HOST data_t data(std::index_sequence<ints...> /*seq*/) {
        return detail::make_tuple<tuple_type>(vecmem::data::vector_view<Ts>(
            vecmem::get_data(detail::get<ints>(_tuple)))...);
    }

    template <typename vec_tuple_data_t, std::size_t... ints>
    DETRAY_DEVICE vtuple::tuple<vector_type<Ts>...> device(
        vec_tuple_data_t& data, std::index_sequence<ints...> /*seq*/) {
        return vtuple::make_tuple(
            vector_type<Ts>(detail::get<ints>(data._tuple))...);
    }
};

/** Tuple of vecmem::vector_view
 */
template <typename vec_tuple_t>
struct vec_tuple_data {
    typename vec_tuple_t::data_t _tuple;

    /** Constructor with vec_tuple - only for host
     */
    template <std::size_t... ints>
    DETRAY_HOST vec_tuple_data(vec_tuple_t& vtuple,
                               std::index_sequence<ints...> seq)
        : _tuple(vtuple.data(seq)) {}
};

/** Get vec_tuple_data
 **/
template <template <typename...> class tuple_type,
          template <typename...> class vector_type, typename... Ts>
inline vec_tuple_data<vec_tuple<tuple_type, vector_type, Ts...> > get_data(
    vec_tuple<tuple_type, vector_type, Ts...>& vtuple) {
    return vec_tuple_data<vec_tuple<tuple_type, vector_type, Ts...> >(
        vtuple, std::make_index_sequence<sizeof...(Ts)>{});
}

// Test function to copy the contenst of vec_tuple_data into vecmem vector
template <template <typename...> class tuple_type>
void tuple_copy(
    vec_tuple_data<vec_tuple<tuple_type, dvector, int, float, double> >& data,
    vecmem::data::vector_view<int>& output1_data,
    vecmem::data::vector_view<float>& output2_data,
    vecmem::data::vector_view<double>& output3_data);

}  // namespace detray
