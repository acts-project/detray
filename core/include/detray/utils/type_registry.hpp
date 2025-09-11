/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

/// @brief match types with indices and vice versa.
///
/// @tparam IDs enum that references the types (not used in base class)
/// @tparam registered_types the types that can be mapped to indices
template <class ID, typename... registered_types>
class type_registry {
    public:
    /// Make the type IDs accessible
    using id = ID;
    /// Make the registered types accessible
    using types = detray::types::list<registered_types...>;

    /// Conventions for some basic info
    enum : std::size_t {
        n_types = sizeof...(registered_types),
        e_any = sizeof...(registered_types),
        e_unknown = sizeof...(registered_types) + 1,
    };

    /// @returns the ID for a type. Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval ID get_id(const object_t& /*obj*/) {
        return get_id<object_t>();
    }

    /// @returns the ID for a type.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval ID get_id() {
        return to_id(detray::types::position<std::decay_t<object_t>, types>);
    }

    /// @returns whether a given types is known in the registry.
    /// Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval bool contains(const object_t& /*obj*/) {
        return contains<object_t>();
    }

    /// @returns whether a given types is known in the registry.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval bool contains() {
        return detray::types::contains<std::decay_t<object_t>, types>;
    }

    /// @returns whether a given index can be mapped to a type.
    DETRAY_HOST_DEVICE static constexpr bool is_valid(
        const std::size_t type_id) {
        return type_id < n_types;
    }

    /// @returns whether a given ID can be mapped to a type.
    template <typename id_t = id>
        requires(!std::convertible_to<id_t, std::size_t>)
    DETRAY_HOST_DEVICE static constexpr bool is_valid(const id_t type_id) {
        return static_cast<std::size_t>(type_id) < n_types;
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
        consteval bool operator()() const noexcept { return is_valid(value); }
    };

    /// Return a type for an index. If the index cannot be mapped, there will be
    /// a compiler error.
    template <ID type_id>
    struct get_type {
        using type = detray::types::at<types, to_index(type_id)>;
    };
};

template <typename registry_t, class type_selector, std::size_t I = 0,
          typename... Fs>
consteval auto map_types(
    const types::list<Fs...>& list = detray::types::list<>{}) {
    using list_t = types::list<Fs...>;

    if constexpr (I == registry_t::n_types) {
        return list;
    } else {
        using next_type = typename registry_t::template get_type<
            static_cast<typename registry_t::id>(I)>::type;
        // Map the mask type to another type
        using mapped_t = typename type_selector::template type<next_type>;

        // Map mask position in mask store to frame type id
        if constexpr (types::contains<mapped_t, list_t>) {
            return map_types<registry_t, type_selector, I + 1u>(list);
        } else {
            return map_types<registry_t, type_selector, I + 1u>(
                types::push_back<list_t, mapped_t>{});
        }
    }
}

template <typename registry_t, class type_selector, std::size_t I = 0,
          typename... Fs>
consteval auto map_ids(const types::list<Fs...>& list,
                       std::array<dindex, registry_t::n_types> id_array = {0}) {
    using list_t = types::list<Fs...>;

    if constexpr (I == registry_t::n_types) {
        return id_array;
    } else {
        using next_type = typename registry_t::template get_type<
            static_cast<typename registry_t::id>(I)>::type;
        // Map the mask type to another type
        using mapped_t = typename type_selector::template type<next_type>;

        // Map mask position in mask store to frame type id
        if constexpr (types::contains<mapped_t, list_t>) {
            id_array[I] = types::position<mapped_t, list_t>;

            return map_ids<registry_t, type_selector, I + 1u>(list, id_array);
        } else {
            id_array[I] = sizeof...(Fs);

            return map_ids<registry_t, type_selector, I + 1u>(
                types::push_back<list_t, mapped_t>{}, id_array);
        }
    }
}

/// @brief Select a number of types from a registry and map them to
/// corresponding IDs
template <typename registry_t, class type_selector_t>
class mapped_type_registry {
    public:
    /// Make the type ids accessible
    using id = typename registry_t::id;
    /// Make the registered types accessible
    using types = decltype(map_types<registry_t, type_selector_t>());

    /// Conventions for some basic info
    enum : std::size_t {
        n_types = detray::types::size<types>,
        e_any = detray::types::size<types>,
        e_unknown = detray::types::size<types> + 1u,
    };

    DETRAY_HOST_DEVICE
    static constexpr const auto& ids() { return id_map; }

    DETRAY_HOST_DEVICE
    static constexpr id mapped_id(id i) {
        return static_cast<id>(id_map[static_cast<std::size_t>(i)]);
    }

    DETRAY_HOST_DEVICE
    static constexpr std::size_t mapped_idx(id i) {
        return id_map[static_cast<std::size_t>(i)];
    }

    /// @returns the ID for a type.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval id get_id() {
        return to_id(detray::types::position<std::decay_t<object_t>, types>);
    }

    /// @returns the ID for a type. Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval id get_id(const object_t&) {
        return get_id<object_t>();
    }

    /// @returns whether a given types is known in the registry.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval bool contains() {
        return detray::types::contains<std::decay_t<object_t>, types>;
    }

    /// @returns whether a given types is known in the registry.
    /// Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval bool contains(const object_t&) {
        return contains<object_t>();
    }

    /// @returns whether a given index can be mapped to a type.
    DETRAY_HOST_DEVICE static constexpr bool is_valid(
        const std::size_t type_idx) {
        if (type_idx >= id_map.size()) {
            return false;
        } else {
            return mapped_idx(static_cast<id>(type_idx)) < n_types;
        }
    }

    /// @returns whether a given ID can be mapped to a type.
    template <typename id_t = id>
        requires(!std::convertible_to<id_t, std::size_t>)
    DETRAY_HOST_DEVICE static constexpr bool is_valid(const id_t type_id) {
        if (static_cast<std::size_t>(type_id) >= id_map.size()) {
            return false;
        } else {
            return mapped_idx(type_id) < n_types;
        }
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr id to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return static_cast<id>(ref_idx);
        }
        if constexpr (ref_idx < detray::types::size<types> - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<id>(detray::types::size<types>);
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr std::size_t to_index(const id i) {
        if (to_id(ref_idx) == i) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return ref_idx;
        }
        if constexpr (ref_idx < detray::types::size<types> - 1) {
            return to_index<ref_idx + 1>(i);
        }
        // This produces a compiler error when used in type unrolling code
        return detray::types::size<types>;
    }

    /// @returns and index for a type
    template <typename object_t>
    struct get_index {
        static constexpr id value = get_id<object_t>();
        DETRAY_HOST_DEVICE
        consteval bool operator()() const noexcept { return is_valid(value); }
    };

    /// @returns a type for an index. If the index cannot be mapped, there will
    /// be a compiler error.
    template <id type_id>
    struct get_type {
        using type = detray::types::at<types, mapped_idx(type_id)>;
    };

    private:
    static constexpr std::array<dindex, registry_t::n_types> id_map =
        map_ids<registry_t, type_selector_t>(types{});
};

}  // namespace detray
