/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/tuple_container.hpp"
#include "detray/definitions/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/** Tuple based container with array elements.
 *
 * @tparam tuple_t is the type of tuple
 * @tparam array_t is the type of array
 * @tparam id_t is the type of indexing integer
 * @tparam Ts... is the types of tuple elements
 */
template <template <typename...> class tuple_t,
          template <typename, std::size_t> class array_t, typename id_t,
          typename... Ts>
class tuple_array_container
    : public tuple_container<tuple_t, id_t,
                             array_t<typename Ts::object_type, Ts::N>...> {

    public:
    // Convenient type declarations
    using base_type =
        tuple_container<tuple_t, id_t,
                        array_t<typename Ts::object_type, Ts::N>...>;
    using base_type::base_type;
    using id_type = typename base_type::id_type;

    static constexpr std::size_t m_tuple_size = base_type::m_tuple_size;

    template <typename... Args>
    using tuple_type = typename base_type::template tuple_type<Args...>;

    template <typename T, std::size_t I>
    using array_type = array_t<T, I>;

    using container_type = typename base_type::container_type;
    using container_data_type =
        tuple_type<array_type<typename Ts::data_type, Ts::N>...>;
    using container_view_type =
        tuple_type<array_type<typename Ts::view_type, Ts::N>...>;

    /**
     * Constructor with a vecmem memory resource. (host-side only)
     *
     * @param resource is the vecmem memory resource
     */
    DETRAY_HOST
    tuple_array_container(vecmem::memory_resource& resource)
        : base_type(container_type(to_host<typename Ts::object_type>(
              resource, std::make_index_sequence<Ts::N>{})...)) {}

    /**
     * Initialize an array with vecmem resource.
     *
     * @tparam T is the type of array elements
     * @tparam ints is the list of index sequence
     *
     * @param resource is the vecmem memory resource
     * @param seq is the index seqeunce (<0,1,...,array_size-1>)
     *
     * @return the array of host objects
     */
    template <typename T, std::size_t... ints>
    DETRAY_HOST array_type<T, sizeof...(ints)> to_host(
        vecmem::memory_resource& resource,
        std::index_sequence<ints...> /*seq*/) {

        array_type<vecmem::memory_resource*, sizeof...(ints)> resources;

        std::fill(resources.begin(), resources.end(), &resource);

        return array_type<T, sizeof...(ints)>{T{resources[ints]}...};
    }

    /**
     * Constructor with a data container.
     * Mainly used in the device side
     *
     * @tparam container_data_t is the type of input data container
     *
     * @param container_data is the data container
     */
    template <typename container_data_t,
              std::enable_if_t<
                  !std::is_base_of_v<vecmem::memory_resource, container_data_t>,
                  bool> = true>
    DETRAY_HOST_DEVICE tuple_array_container(container_data_t& container_data)
        : base_type(initialize_device_arrays(
              container_data, std::make_index_sequence<m_tuple_size>{})) {}

    /**
     * Generate the tuple of device array
     *
     * @tparam container_data_t is the type of input data container
     * @tparam ints is the list of index sequence
     *
     * @param container_data is the data container
     * @param seq is the index seqeunce (<0,1,...,tuple_size-1>)
     *
     * @return the tuple of device array
     */
    template <typename container_data_t, std::size_t... ints>
    DETRAY_HOST_DEVICE container_type
    initialize_device_arrays(container_data_t& container_data,
                             std::index_sequence<ints...> /*seq*/) {

        return container_type(to_device<Ts::object_type>(
            detail::get<ints>(container_data.m_data))...);
    }

    /**
     * Convert the data objects of array into device objects
     *
     * @tparam object_t is the type of device object
     * @tparam data_t is the type of data object
     * @tparam N is the size of array
     *
     * @param A is the input data array
     *
     * @return the array of device objects
     */
    template <typename object_t, typename data_t, std::size_t N>
    DETRAY_HOST_DEVICE array_type<object_t, N> to_device(
        array_type<data_t, N>& A) {
        return to_device<object_t>(A, std::make_index_sequence<N>{});
    }

    /**
     * Convert the data objects of array into device objects
     *
     * @tparam object_t is the type of device object
     * @tparam data_t is the type of data object
     * @tparam N is the size of array
     * @tparam ints is the list of index sequence
     *
     * @param A is the input data array
     * @param seq is the index seqeunce (<0,1,...,array_size-1>)
     *
     * @return the array of device objects
     */
    template <typename object_t, typename data_t, std::size_t N,
              std::size_t... ints>
    DETRAY_HOST_DEVICE array_type<object_t, N> to_device(
        array_type<data_t, N>& A, std::index_sequence<ints...> /*seq*/) {
        return array_type<object_t, N>{object_t(A[ints])...};
    }

    /**
     * Generate the tuple of device array
     *
     * @tparam container_data_t is the type of input data container
     * @tparam ints is the list of index sequence
     *
     * @param container_data is the data container
     * @param seq is the index seqeunce (<0,1,...,tuple_size-1>)
     *
     * @return the tuple of device array
     */
    template <std::size_t... ints>
    DETRAY_HOST container_data_type
    initialize_data_arrays(vecmem::memory_resource& resource,
                           std::index_sequence<ints...> /*seq*/) {

        return container_data_type{(to_data<typename Ts::data_type>(
            detail::get<ints>(this->m_container), resource))...};
    }

    /**
     * Convert the host objects of array into data objects
     *
     * @tparam object_t is the type of host object
     * @tparam data_t is the type of data object
     * @tparam N is the size of array
     *
     * @param A is the input host array
     *
     * @return the array of data objects
     */
    template <typename data_t, typename object_t, std::size_t N>
    DETRAY_HOST array_type<data_t, N> to_data(
        array_type<object_t, N>& A, vecmem::memory_resource& resource) {
        return to_data<data_t>(A, resource, std::make_index_sequence<N>{});
    }

    /**
     * Convert the host objects of array into data objects
     *
     * @tparam object_t is the type of host object
     * @tparam data_t is the type of data object
     * @tparam N is the size of array
     * @tparam ints is the list of index sequence
     *
     * @param A is the input host array
     * @param seq is the index seqeunce (<0,1,...,array_size-1>)
     *
     * @return the array of data objects
     */
    template <typename data_t, typename object_t, std::size_t N,
              std::size_t... ints>
    DETRAY_HOST array_type<data_t, N> to_data(
        array_type<object_t, N>& A, vecmem::memory_resource& resource,
        std::index_sequence<ints...> /*seq*/) {

        (void)resource;
        return array_type<data_t, N>{vecmem::get_data(A[ints])...};
    }
};

/**
 * A static implementation of data container for device
 *
 * @tparam container_t is the type for host container
 */
template <typename container_t>
struct tuple_array_container_data {

    template <typename... Args>
    using tuple_type = typename container_t::template tuple_type<Args...>;

    using container_view_type = typename container_t::container_view_type;

    /**
     * Constructor with a tuple array container (host-side only)
     *
     * @param container is the input container
     * @param resource is the vecmem memory resource
     */
    DETRAY_HOST
    tuple_array_container_data(container_t& container,
                               vecmem::memory_resource& resource)
        : m_data(container.initialize_data_arrays(
              resource, std::make_index_sequence<container.size()>{})) {}

    typename container_t::container_data_type m_data;
};

/**
 * A static implementation of view container for device
 *
 * @tparam container_t is the type for host container
 */
template <typename container_t>
struct tuple_array_container_view {

    using container_view_type = typename container_t::container_view_type;

    /**
     * Constructor with a tuple array data container (host-side only)
     *
     * @param container_data is the input data container
     */
    DETRAY_HOST
    tuple_array_container_view(
        tuple_array_container_data<container_t>& container_data)
        : m_data(initialize_view_arrays(
              container_data,
              std::make_index_sequence<container_t::m_tuple_size>{})) {}

    /**
     * Generate the tuple of view array
     *
     * @tparam ints is the list of index sequence
     *
     * @param seq is the index seqeunce (<0,1,...,tuple_size-1>)
     *
     * @return the tuple of view array
     */
    template <std::size_t... ints>
    DETRAY_HOST container_view_type initialize_view_arrays(
        tuple_array_container_data<container_t>& container_data,
        std::index_sequence<ints...> /*seq*/) {
        return container_view_type(detail::get<ints>(container_data.m_data)...);
    }

    container_view_type m_data;
};

/**
 * A stand-alone function to get the data container
 *
 * @tparam container_t is the type for host container
 *
 * @return the tuple array data container
 */
template <typename container_t>
inline tuple_array_container_data<container_t> get_data(
    container_t& container, vecmem::memory_resource& resource) {
    return {container, resource};
}

}  // namespace detray