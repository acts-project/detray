/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <iostream>
#include <type_traits>

#include "detray/core/mask_store.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Base class for a type registry that allows to map indices to types and vice
 * versa.
 *
 * @tparam registered_types the types that can be mapped to indices
 */
template <typename... registered_types>
class registry_base {
    public:
    /** Conventions for some basic info */
    enum : unsigned int {
        n_types = sizeof...(registered_types),
        e_any = 0,  // sizeof...(registered_types),
        e_unknown = sizeof...(registered_types) + 1,
    };

    /** Get the index for a type. Needs to be unrolled in case of thrust tuple.
     */
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr unsigned int get_id() {
        return unroll_ids<object_t, registered_types...>();
    }

    /** Checks whether a given types is known in the registry.*/
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr bool is_defined() {
        return not(get_id<object_t>() == e_unknown);
    }

    /** Checks whether a given index can be mapped to a type.*/
    template <unsigned int ID>
    DETRAY_HOST_DEVICE static constexpr bool is_valid_id() {
        return 0 <= ID and ID < n_types;
    }

    /** Extract an index and check it.*/
    template <typename object_t>
    struct get_index {
        static constexpr unsigned int value = get_id<object_t>();
        constexpr bool operator()() noexcept { return is_valid_id<value>(); }
    };

    /** Return a type for an index. If the index cannot be mapped, there will be
     * a compiler error.
     */
    template <unsigned int index, template <typename...> class tuple_t = dtuple>
    struct get_type {
        using type = decltype(std::get<index>(tuple_t<registered_types...>{}));
    };

    private:
    /// dummy type
    struct empty_type {};

    /** Gets the position of a type in a parameter pack, without using tuples.*/
    template <typename object_t, typename first_t = empty_type,
              typename... remaining_types>
    DETRAY_HOST_DEVICE static constexpr unsigned int unroll_ids() {
        if constexpr (not std::is_same_v<first_t, empty_type> and
                      not std::is_same_v<object_t, first_t>) {
            return unroll_ids<object_t, remaining_types...>();
        }
        if constexpr (std::is_same_v<object_t, first_t>) {
            return n_types - sizeof...(remaining_types) - 1;
        }
        return e_unknown;
    }
};

/** Registry for geometric objects.*/
template <typename... registered_types>
class default_object_registry : public registry_base<registered_types...> {
    public:
    using type_registry = registry_base<registered_types...>;

    enum : unsigned int {
        n_types = type_registry::n_types,
    };

    /// Known primitives
    enum id : unsigned int {
        e_surface = 0,
        e_portal = 0,  // not used (same as surface)
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <unsigned int ID, template <typename...> class tuple_t = dtuple>
    using get_type = typename type_registry::template get_type<ID, tuple_t>;

    /*template <unsigned int current_id = 0>
    DETRAY_HOST_DEVICE static auto &get_dynamic_type(const unsigned int index) {
        if (current_id == index) {
            get_type<current_id> type_wrapper{};
            return type_wrapper;
        }
        // Next type
        if constexpr (current_id < n_types - 1) {
            return get_dynamic_type<current_id + 1>(index);
        }
    }*/
};

/** Registry object for masks */
template <typename... registered_types>
class default_mask_registry : public registry_base<registered_types...> {
    public:
    using type_registry = registry_base<registered_types...>;

    enum : unsigned int {
        n_types = type_registry::n_types,
    };

    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_container = mask_store<tuple_t, vector_t, registered_types...>;
    using mask_link = typename mask_container<>::mask_link;
    using mask_range = typename mask_container<>::mask_range;

    /// Give your mask types a name (needs to be consecutive to be matched)
    /// to a type!)
    enum id : unsigned int {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_portal_cylinder3 = 3,  // no distinction from surface cylinder
        e_ring2 = 4,
        e_portal_ring2 = 4,
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <unsigned int ID, template <typename...> class tuple_t = dtuple>
    using get_type = typename type_registry::template get_type<ID, tuple_t>;

    /*template <unsigned int current_id = 0>
    DETRAY_HOST_DEVICE static auto &get_dynamic_type(const unsigned int index) {
        if (current_id == index) {
            return get_type<current_id>{};
        }
        // Next mask type
        if constexpr (current_id < mask_container<>::size()) {
            return get_dynamic_type<current_id + 1>(index);
        }
    }*/

    /*template <unsigned int current_type_index = e_surface, typename
    args_wrapper_t = empty_type, typename first_t = empty_type, typename
    ...remaining_types> auto type_matcher(const unsigned int index,
    args_wrapper_t& args, Args &&... args) { if constexpr (not
    std::is_same_v<first_t, empty_type> and is_defined<current_type_index>()) {
            if (index != current_type_index) {
                type_matcher<current_type_index + 1>(index, args, ...);
            }
        }
        using callable = typename get_type<current_type_index>::type;
        return callable(args...);
    }*/
};

/** Registry class for surface finders (e.g. grids) */
template <typename... registered_types>
class default_sf_finder_registry : public registry_base<registered_types...> {
    public:
    using type_registry = registry_base<registered_types...>;

    enum : unsigned int {
        n_types = type_registry::n_types,
    };

    /// Surface finders
    enum id : unsigned int {
        e_circ_reg_grid = 0,
        e_reg_reg_grid = 1,  // not used (same as surface)
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <unsigned int ID, template <typename...> class tuple_t = dtuple>
    using get_type = typename type_registry::template get_type<ID, tuple_t>;

    /*DETRAY_HOST_DEVICE static constexpr auto get_tp(const unsigned int index)
    { return unroll_types<registered_types...>(index);
    }

    template<template <typename...> class tuple_t = dtuple>
    struct match_type {
        unsigned int _index;
        match_type(unsigned int index) : _index(index) {}
        using type = decltype(get_tp(_index));
    };

    template <typename first_t = empty_type, typename ...remaining_types>
    DETRAY_HOST_DEVICE static auto unroll_types (const unsigned int index) ->
    get_type<n_types - sizeof...(remaining_types) - 1> { constexpr unsigned int
    current_type_index = n_types - sizeof...(remaining_types) - 1; if (not
    std::is_same_v<first_t, empty_type> and not (index == current_type_index)) {
            unroll_types<remaining_types...>(index);
        }
        if (index == current_type_index) {
            return get_type<current_type_index>{};
        }
        return get_type<e_unknown>{};
    }*/
};

}  // namespace detray