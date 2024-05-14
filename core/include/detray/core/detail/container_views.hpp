/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/tuple_helpers.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/// Container types used in device code
using device_container_types =
    container_types<vecmem::device_vector, detray::tuple, darray,
                    vecmem::jagged_device_vector>;

namespace detail {

// Views for types that aggregate containers/other viewable types

/// Empty view type for inheritance template resolution
struct dbase_view {};

/// Container view helper that aggregates multiple vecmem views and performs
/// compile-time checks.
template <bool /*check value*/, typename... view_ts>
class dmulti_view_helper {};

/// In case the checks fail
template <typename... view_ts>
class dmulti_view_helper<false, view_ts...> {};

/// @brief General view type that aggregates vecmem based view implementations.
///
/// This is for detray types that hold multiple members that all define custom
/// view types of their own. The 'sub'-views are being aggregated in this helper
/// and are extracted in the types contructor and then handed down to the
/// member constructors.
template <typename... view_ts>
struct dmulti_view_helper<true, view_ts...> : public dbase_view {
    dtuple<view_ts...> m_view;

    dmulti_view_helper() = default;

    /// Tie multiple views together
    DETRAY_HOST_DEVICE
    explicit dmulti_view_helper(view_ts&&... views)
        : m_view(std::move(views)...) {}

    /// Tie multiple views together
    DETRAY_HOST_DEVICE
    explicit dmulti_view_helper(view_ts&... views) : m_view(views...) {}
};

/// Helper trait to determine if a type can be interpreted as a (composite)
/// vecemem view
/// @{
template <typename T, typename = void>
struct is_device_view : public std::false_type {};

template <typename T>
struct is_device_view<
    T, std::enable_if_t<
           std::is_base_of_v<detray::detail::dbase_view,
                             std::remove_reference_t<std::remove_cv_t<T>>>,
           void>> : public std::true_type {};

template <typename T>
inline constexpr bool is_device_view_v = is_device_view<T>::value;
/// @}

/// Helper trait to check whether a type has a vecmem view defined
/// @{
template <class T, typename = void>
struct has_view : public std::false_type {
    using type = void;
};

template <class T>
struct has_view<
    T, std::enable_if_t<detray::detail::is_device_view_v<typename T::view_type>,
                        void>> : public std::true_type {
    using type =
        std::conditional_t<std::is_const_v<T>, typename T::const_view_type,
                           typename T::view_type>;
};

template <class T>
inline constexpr bool has_view_v = has_view<T>::value;

template <class T>
using get_view_t = typename has_view<T>::type;
/// @}

}  // namespace detail

/// The detray container view exists, if all contained view types also derive
/// from @c dbase_view.
template <typename... view_ts>
using dmulti_view = detray::detail::dmulti_view_helper<
    std::conjunction_v<detail::is_device_view<view_ts>...>, view_ts...>;

/// @brief Detray version of 'get_data' - non-const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <class T,
          std::enable_if_t<detail::is_device_view_v<typename T::view_type>,
                           bool> = true>
typename T::view_type get_data(T& viewable) {
    return viewable.get_data();
}

/// @brief Detray version of 'get_data' - const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <class T,
          std::enable_if_t<detail::is_device_view_v<typename T::view_type>,
                           bool> = true>
typename T::const_view_type get_data(const T& viewable) {
    return viewable.get_data();
}

/// Type trait specializations for vecmem containers
/// @{

/// Specialized view for @c vecmem::vector containers
template <typename T>
using dvector_view = vecmem::data::vector_view<T>;

/// Get the vecmem view type of a vector
template <typename T, typename A>
dvector_view<T> get_data(std::vector<T, A>& vec) {
    return vecmem::get_data(vec);
}

/// Get the vecmem view type of a vector - const
template <typename T, typename A>
dvector_view<const T> get_data(const std::vector<T, A>& vec) {
    return vecmem::get_data(vec);
}

/// Specialized view for @c vecmem::jagged_vector containers
template <typename T>
using djagged_vector_view = vecmem::data::jagged_vector_view<T>;

/// Specialization of 'is view' for @c vecmem::data::vector_view containers
template <typename T>
struct detail::is_device_view<vecmem::data::vector_view<T>, void>
    : public std::true_type {};

/// Specialization of 'is view' for constant @c vecmem::data::vector_view
/// containers
template <typename T>
struct detail::is_device_view<const vecmem::data::vector_view<T>, void>
    : public std::true_type {};

/// Specialization of the view getter for @c vecmem::vector
template <typename T>
struct detail::has_view<vecmem::vector<T>, void> : public std::true_type {
    using type = dvector_view<T>;
};

/// Specialization of the view getter for @c vecmem::vector - const
template <typename T>
struct detail::has_view<const vecmem::vector<T>, void> : public std::true_type {
    using type = dvector_view<const T>;
};

/// Specialization of the view getter for @c vecmem::device_vector
template <typename T>
struct detail::has_view<vecmem::device_vector<T>, void>
    : public std::true_type {
    using type = dvector_view<T>;
};

/// Specialization of the view getter for @c vecmem::device_vector - const
template <typename T>
struct detail::has_view<const vecmem::device_vector<T>, void>
    : public std::true_type {
    using type = dvector_view<const T>;
};

/// Specialized view for @c vecmem::jagged_vector containers
template <typename T>
using djagged_vector_view = vecmem::data::jagged_vector_view<T>;

/// Specialization of 'is view' for @c vecmem::data::jagged_vector_view
/// containers
template <typename T>
struct detail::is_device_view<vecmem::data::jagged_vector_view<T>, void>
    : public std::true_type {};

/// Specialization of 'is view' for constant @c vecmem::data::jagged_vector_view
/// containers
template <typename T>
struct detail::is_device_view<const vecmem::data::jagged_vector_view<T>, void>
    : public std::true_type {};

/// Specialization of the view getter for @c vecmem::jagged_vector
template <typename T>
struct detail::has_view<vecmem::jagged_vector<T>, void>
    : public std::true_type {
    using type = djagged_vector_view<T>;
};

/// Specialization of the view getter for @c vecmem::jagged_vector - const
template <typename T>
struct detail::has_view<const vecmem::jagged_vector<T>, void>
    : public std::true_type {
    using type = djagged_vector_view<const T>;
};

/// Specialization of the view getter for @c vecmem::jagged_device_vector
template <typename T>
struct detail::has_view<vecmem::jagged_device_vector<T>, void>
    : public std::true_type {
    using type = djagged_vector_view<T>;
};

/// Specialization of the view getter for @c vecmem::jagged_device_vector
/// - const
template <typename T>
struct detail::has_view<const vecmem::jagged_device_vector<T>, void>
    : public std::true_type {
    using type = djagged_vector_view<const T>;
};
/// @}

}  // namespace detray
