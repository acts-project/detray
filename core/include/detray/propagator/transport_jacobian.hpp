/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "./codegen/transport_jacobian.hpp"

namespace detray::detail {
template <typename T>
concept is_transport_jacobian_type =
    requires { typename T::algebra_type; } &&
    (std::same_as<T, transport_jacobian_matrix_with_gradient<
                         typename T::algebra_type>> ||
     std::same_as<T, transport_jacobian_matrix_without_gradient<
                         typename T::algebra_type>>);
}  // namespace detray::detail

namespace detray::getter {}

namespace detray::matrix {
template <typename matrix_type>
DETRAY_HOST_DEVICE constexpr auto identity()
    requires detail::is_transport_jacobian_type<matrix_type>
{
    return matrix_type::identity();
}
}  // namespace detray::matrix

namespace algebra::traits {
template <typename algebra_t>
struct dimensions<
    ::detray::detail::transport_jacobian_matrix_with_gradient<algebra_t>> {

    using size_type = std::size_t;

    static constexpr size_type dim{2};
    static constexpr size_type rows{8};
    static constexpr size_type columns{8};
};

template <typename algebra_t>
struct dimensions<
    ::detray::detail::transport_jacobian_matrix_without_gradient<algebra_t>> {

    using size_type = std::size_t;

    static constexpr size_type dim{2};
    static constexpr size_type rows{8};
    static constexpr size_type columns{8};
};

template <typename algebra_t>
struct index<
    ::detray::detail::transport_jacobian_matrix_with_gradient<algebra_t>> {

    using type = std::size_t;
};

template <typename algebra_t>
struct index<
    ::detray::detail::transport_jacobian_matrix_without_gradient<algebra_t>> {

    using type = std::size_t;
};
}  // namespace algebra::traits
