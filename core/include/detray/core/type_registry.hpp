/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>
#include <utility>

#include "detray/core/mask_store.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

/** Base class for a type registry that allows to map indices to types and vice
 * versa.
 *
 * @tparam IDs enum that references the types (not used in base class)
 * @tparam registered_types the types that can be mapped to indices
 */
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
    /** Conventions for some basic info */
    enum : unsigned int {
        n_types = sizeof...(registered_types),
        e_any = sizeof...(registered_types),
        e_unknown = sizeof...(registered_types) + 1,
    };

    /** Get the index for a type. Needs to be unrolled in case of thrust tuple.
     */
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr unsigned int get_id() {
        return unroll_ids<std::remove_reference_t<object_t>,
                          registered_types...>();
    }

    /** Checks whether a given types is known in the registry.*/
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr bool is_defined() {
        return not(get_id<object_t>() == e_unknown);
    }

    /** Checks whether a given index can be mapped to a type.*/
    template <ID type_id>
    DETRAY_HOST_DEVICE static constexpr bool is_valid_id() {
        return 0 <= type_id and type_id < n_types;
    }

    /** Extract an index and check it.*/
    template <typename object_t>
    struct get_index {
        static constexpr ID value = get_id<object_t>();
        constexpr bool operator()() noexcept { return is_valid_id<value>(); }
    };

    /** Return a type for an index. If the index cannot be mapped, there will be
     * a compiler error.
     */
    template <ID type_id, template <typename...> class tuple_t = dtuple>
    struct get_type {
        using type = std::remove_reference_t<decltype(std::get<type_id>(
            tuple_t<registered_types...>{}))>;
    };

    private:
    /// dummy type
    struct empty_type {};

    /** Gets the position of a type in a parameter pack, without using tuples.*/
    template <typename object_t, typename first_t = empty_type,
              typename... remaining_types>
    DETRAY_HOST_DEVICE static constexpr ID unroll_ids() {
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
class object_registry
    : public registry_base<unsigned int, true, registered_types...> {
    public:
    using type_registry =
        registry_base<unsigned int, true, registered_types...>;

    enum : unsigned int {
        n_types = type_registry::n_types,
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    /// Known primitives
    enum id : unsigned int {
        e_surface = 0,
        e_portal = 0,  // not used (same as surface)
    };

    using link_type = dindex;
    using range_type = dindex_range;

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <unsigned int type_id,
              template <typename...> class tuple_t = dtuple>
    using get_type =
        typename type_registry::template get_type<type_id, tuple_t>;
};

/** Registry object for masks */
template <class ID, typename... registered_types>
class mask_registry
    : public registry_base<
          ID, std::is_enum_v<ID> and std::is_convertible_v<ID, unsigned int>,
          registered_types...> {
    public:
    using type_registry = registry_base<
        ID, std::is_enum_v<ID> and std::is_convertible_v<ID, unsigned int>,
        registered_types...>;

    enum : unsigned int {
        n_types = type_registry::n_types,
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    // Make the type IDs accessible
    using id = ID;

    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using container_type = mask_store<tuple_t, vector_t, registered_types...>;
    // using link_type = typed_index<id>;
    using link_type = typename container_type<>::link_type;
    using range_type = typename container_type<>::range_type;

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <ID type_id, template <typename...> class tuple_t = dtuple>
    using get_type =
        typename type_registry::template get_type<type_id, tuple_t>;
};

/** Registry class for surface finders (e.g. grids) */
// TODO: Merge with mask registry
template <typename... registered_types>
class sf_finder_registry
    : public registry_base<unsigned int, true, registered_types...> {
    public:
    using type_registry =
        registry_base<unsigned int, true, registered_types...>;

    enum : unsigned int {
        n_types = type_registry::n_types,
        e_any = type_registry::e_any,
        e_unknown = type_registry::e_unknown,
    };

    /// Surface finders
    enum id : std::size_t {
        e_brute_force = 0,
        e_z_phi_grid = 1,  // barrel
        e_r_phi_grid = 2,  // endcap
    };

    // using link_type = typed_index<id>;
    // using range_type = typed_index<id, dindex_range>;
    using link_type = std::array<dindex, 2>;
    using range_type = dindex_range;

    template <typename T>
    using get_index = typename type_registry::template get_index<T>;

    template <unsigned int type_id,
              template <typename...> class tuple_t = dtuple>
    using get_type =
        typename type_registry::template get_type<type_id, tuple_t>;
};

}  // namespace detray