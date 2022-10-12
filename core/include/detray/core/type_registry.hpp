/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/tuple_array_container.hpp"
#include "detray/core/detail/tuple_vector_container.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

/// Base class for a type registry that allows to map indices to types and vice
/// versa.
///
/// @tparam IDs enum that references the types (not used in base class)
/// @tparam registered_types the types that can be mapped to indices
// TODO: Merge with tuple_container
template <class ID, bool /*put checks on IDs type*/,
          typename... registered_types>
class registry_base;

template <class ID, typename... registered_types>
class registry_base<ID, false, registered_types...> {
    // Produce meaningful errors
    static_assert(
        std::is_enum_v<ID>,
        "First template parameter of a type registry must be the type enum!");
    static_assert(std::is_convertible_v<ID, std::size_t>,
                  "Type enum must be convertible to std::size_t!");
};

template <class ID, typename... registered_types>
class registry_base<ID, true, registered_types...> {
    public:
    /// Conventions for some basic info
    enum : std::size_t {
        n_types = sizeof...(registered_types),
        e_any = sizeof...(registered_types),
        e_unknown = sizeof...(registered_types) + 1,
    };

    /// Get the index for a type. Needs to be unrolled in case of thrust tuple.
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr ID get_id() {
        return unroll_ids<std::remove_reference_t<object_t>,
                          registered_types...>();
    }

    /// Get the index for a type. Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr ID get_id(object_t& /*obj*/) {
        return get_id<object_t>();
    }

    /// Checks whether a given types is known in the registry.
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr bool is_defined() {
        return not(get_id<object_t>() == e_unknown);
    }

    /// Checks whether a given index can be mapped to a type.
    DETRAY_HOST_DEVICE static constexpr bool is_valid(
        const std::size_t type_id) {
        return type_id < n_types;
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr ID to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return static_cast<ID>(ref_idx);
        }
        if constexpr (ref_idx < sizeof...(registered_types) - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<ID>(sizeof...(registered_types));
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr std::size_t to_index(const ID id) {
        if (to_id(ref_idx) == id) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return ref_idx;
        }
        if constexpr (ref_idx < sizeof...(registered_types) - 1) {
            return to_index<ref_idx + 1>(id);
        }
        // This produces a compiler error when used in type unrolling code
        return sizeof...(registered_types);
    }

    /// Extract an index and check it.
    template <typename object_t>
    struct get_index {
        static constexpr ID value = get_id<object_t>();
        DETRAY_HOST_DEVICE
        constexpr bool operator()() noexcept { return is_valid(value); }
    };

    /// Return a type for an index. If the index cannot be mapped, there will be
    /// a compiler error.
    template <ID type_id, template <typename...> class tuple_t = dtuple>
    struct get_type {
        using type = std::remove_reference_t<decltype(
            detail::get<to_index(type_id)>(tuple_t<registered_types...>{}))>;
    };

    private:
    /// dummy type
    struct empty_type {};

    /// Gets the position of a type in a parameter pack, without using tuples.
    template <typename object_t, typename first_t = empty_type,
              typename... remaining_types>
    DETRAY_HOST_DEVICE static constexpr ID unroll_ids() {
        if constexpr (not std::is_same_v<first_t, empty_type> and
                      not std::is_same_v<object_t, first_t>) {
            return unroll_ids<object_t, remaining_types...>();
        }
        if constexpr (std::is_same_v<object_t, first_t>) {
            return static_cast<ID>(n_types - sizeof...(remaining_types) - 1);
        }
        return static_cast<ID>(e_unknown);
    }
};

/// Tuple vector container registry
template <class ID, typename... registered_types>
class tuple_vector_registry
    : public registry_base<ID, std::is_enum_v<ID>, registered_types...> {
    public:
    using type_registry =
        registry_base<ID, std::is_enum_v<ID>, registered_types...>;

    enum : std::size_t {
        n_types = type_registry::n_types,
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    // Make the type IDs accessible
    using id = ID;

    // Cuda cannot handle ID non-types here, so leave it for now
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using store_type =
        tuple_vector_container<tuple_t, vector_t, ID, registered_types...>;
    using link_type = dtyped_index<ID, dindex>;
    using range_type = dtyped_index<ID, dindex_range>;

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <ID type_id, template <typename...> class tuple_t = dtuple>
    using get_type =
        typename type_registry::template get_type<type_id, tuple_t>;
};

/// Registry class for surface finders (e.g. grids)
// TODO: Merge with mask registry
template <class ID, typename...>
class tuple_array_registry;

/// Specialization to resolve template parameter packs
template <class ID, std::size_t... sizes, typename... registered_types>
class tuple_array_registry<ID, std::index_sequence<sizes...>,
                           registered_types...>
    : public registry_base<ID, std::is_enum_v<ID>, registered_types...> {
    public:
    using type_registry =
        registry_base<ID, std::is_enum_v<ID>, registered_types...>;

    enum : std::size_t {
        n_types = type_registry::n_types,
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    // Make the type IDs accessible
    using id = ID;

    // Cuda cannot handle ID non-types here, so leave it for now
    template <template <typename...> class tuple_t = dtuple,
              template <typename, std::size_t> class array_t = darray>
    using store_type = tuple_array_container<tuple_t, array_t, ID,
                                             std::index_sequence<sizes...>,
                                             registered_types...>;
    using link_type = dtyped_index<ID, dindex>;
    using range_type = dtyped_index<ID, dindex_range>;

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <ID type_id, template <typename...> class tuple_t = dtuple>
    using get_type =
        typename type_registry::template get_type<type_id, tuple_t>;
};

}  // namespace detray