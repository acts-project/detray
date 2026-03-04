/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray::getter {
namespace concepts {
template <typename T>
concept permits_compile_time_matrix_index_element = requires(const T& t) {
    { t.template element<0, 0>() };
};

template <typename T>
concept permits_compile_time_vector_index_element = requires(const T& t) {
    { t.template element<0>() };
};
}  // namespace concepts

template <std::size_t I, std::size_t J, typename T>
DETRAY_HOST_DEVICE decltype(auto) element(T& matrix) {
    if constexpr (concepts::permits_compile_time_matrix_index_element<T>) {
        return matrix.template element<I, J>();
    } else {
        return element(matrix, I, J);
    }
}

template <std::size_t I, typename T>
DETRAY_HOST_DEVICE decltype(auto) element(T& vector) {
    if constexpr (concepts::permits_compile_time_vector_index_element<T>) {
        return vector.template element<I>();
    } else {
        return element(vector, I);
    }
}
}  // namespace detray::getter
