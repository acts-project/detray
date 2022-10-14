/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/tuple_container.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/** Tuple based container with vector elements.
 *
 * @tparam tuple_t is the type of tuple
 * @tparam vector_t is the type of vector
 * @tparam id_t is an enum that is compared to an indexing integer
 * @tparam Ts are the types of tuple elements
 */
template <template <typename...> class tuple_t,
          template <typename...> class vector_t, typename id_t, typename... Ts>
class tuple_vector_container final
    : public tuple_container<tuple_t, id_t, vector_t<Ts>...> {

    public:
    // Convenient type declarations
    using base_type = tuple_container<tuple_t, id_t, vector_t<Ts>...>;
    using base_type::base_type;
    using id_type = typename base_type::id_type;
    using base_type::to_index;

    static constexpr std::size_t m_tuple_size = base_type::m_tuple_size;

    template <typename... Args>
    using tuple_type = typename base_type::template tuple_type<Args...>;

    template <typename T>
    using vector_type = vector_t<T>;

    using container_type = typename base_type::container_type;
    using container_data_type = tuple_type<vecmem::data::vector_view<Ts>...>;

    /**
     * Constructor with a vecmem memory resource. (host-side only)
     *
     * @param resource is the vecmem memory resource
     */
    DETRAY_HOST
    tuple_vector_container(vecmem::memory_resource &resource)
        : base_type(container_type(vector_type<Ts>{&resource}...)) {}

    /**
     * Constructor with a data container.
     * Mainly used in the device side
     *
     * @tparam is the type of input data container
     *
     * @param container_data is the data container
     */
    template <
        typename container_data_t,
        std::enable_if_t<
            !std::is_base_of_v<vecmem::memory_resource, container_data_t> &&
                !std::is_same_v<tuple_vector_container, container_data_t>,
            bool> = true>
    DETRAY_HOST_DEVICE tuple_vector_container(container_data_t &container_data)
        : base_type(initialize_device_vectors(
              container_data, std::make_index_sequence<m_tuple_size>{})) {}

    /**
     * Generate the tuple of device vector
     *
     * @tparam container_data_t is the type of input data container
     * @tparam ints is the list of index sequence
     *
     * @param container_data is the data container
     * @param seq is the index seqeunce (<0,1,...,tuple_size-1>)
     *
     * @return the tuple of device vector
     */
    template <typename container_data_t, std::size_t... ints>
    DETRAY_HOST_DEVICE container_type
    initialize_device_vectors(container_data_t &container_data,
                              std::index_sequence<ints...> /*seq*/) {
        return vtuple::make_tuple(
            vector_type<Ts>(detail::get<ints>(container_data.m_data))...);
    }

    /** Add a new value in place
     *
     * @tparam ID is the index of target vector
     * @tparam Args are the types of the constructor arguments
     *
     * @param args is the list of constructor argument
     *
     * @note in general can throw an exception
     */
    template <id_t ID, typename... Args>
    DETRAY_HOST auto &add_value(Args &&... args) noexcept(false) {

        auto &gr = detail::get<to_index(ID)>(this->m_container);

        return gr.emplace_back(std::forward<Args>(args)...);
    }

    /** Add a new vector
     *
     * @tparam current_id is the index of target vector
     * @tparam T is the value type
     *
     * @param vec is the vector to be added
     *
     * @note in general can throw an exception
     */
    template <typename T>
    DETRAY_HOST inline void add_vector(vector_type<T> &vec) noexcept(false) {

        static_assert((std::is_same_v<T, Ts> || ...) == true,
                      "The type is not included in the parameter pack.");

        auto &gr = detail::get<vector_type<T>>(this->m_container);

        gr.reserve(gr.size() + vec.size());
        gr.insert(gr.end(), vec.begin(), vec.end());
    }

    /** Add a new vector (move semantics)
     *
     * @tparam current_id is the index of target vector
     * @tparam T is the value type
     *
     * @param vec is the vector to be added
     *
     * @note in general can throw an exception
     */
    template <typename T>
    DETRAY_HOST inline void add_vector(vector_type<T> &&vec) noexcept(false) {

        static_assert((std::is_same_v<T, Ts> || ...) == true,
                      "The type is not included in the parameter pack.");

        auto &gr = detail::get<vector_type<T>>(this->m_container);

        gr.reserve(gr.size() + vec.size());
        gr.insert(gr.end(), std::make_move_iterator(vec.begin()),
                  std::make_move_iterator(vec.end()));
    }

    /** Append a container to the current one
     *
     * @tparam current_idx is the index to start unrolling
     *
     * @param other The other container
     *
     * @note in general can throw an exception
     */
    template <std::size_t current_idx = 0>
    DETRAY_HOST inline void append_container(
        tuple_vector_container &other) noexcept(false) {
        auto &gr = detail::get<current_idx>(other);
        add_vector(gr);

        if constexpr (current_idx < sizeof...(Ts) - 1) {
            append_container<current_idx + 1>(other);
        }
    }
};

/**
 * A static implementation of data container for device
 *
 * @tparam container_t is the type for host container
 */
template <typename container_t>
struct tuple_vector_container_data {

    template <typename... Args>
    using tuple_type = typename container_t::template tuple_type<Args...>;
    using container_type = typename container_t::container_type;
    using container_data_type = typename container_t::container_data_type;

    /**
     * Constructor with a tuple vector container (host-side only)
     *
     * @tparam ints is the list of index sequence
     *
     * @param container is the input tuple vector container
     * @param seq is the index seqeunce (<0,1,...,tuple_size-1>)
     */
    template <std::size_t... ints>
    DETRAY_HOST tuple_vector_container_data(container_t &container,
                                            std::index_sequence<ints...> seq)
        : m_data(initialize_data_vectors(container, seq)) {}

    /**
     * Generate the tuple of data vector
     *
     * @tparam ints is the list of index sequence
     *
     * @param container is the input tuple vector container
     * @param seq is the index seqeunce (<0,1,...,tuple_size-1>)
     *
     * @return the tuple of data vector
     */
    template <std::size_t... ints>
    DETRAY_HOST container_data_type initialize_data_vectors(
        container_t &container, std::index_sequence<ints...> /*seq*/) {

        return detail::make_tuple<tuple_type>(
            vecmem::get_data(detail::get<ints>(container.get()))...);
    }

    container_data_type m_data;
};

/**
 * A stand-alone function to get the data container
 *
 * @tparam container_t is the type for host container
 *
 * @return the tuple vector data container
 */
template <template <typename...> class tuple_t,
          template <typename...> class vector_t, typename id_t, typename... Ts>
inline tuple_vector_container_data<
    tuple_vector_container<tuple_t, vector_t, id_t, Ts...>>
get_data(tuple_vector_container<tuple_t, vector_t, id_t, Ts...> &container) {
    return tuple_vector_container_data<
        tuple_vector_container<tuple_t, vector_t, id_t, Ts...>>(
        container, std::make_index_sequence<tuple_vector_container<
                       tuple_t, vector_t, id_t, Ts...>::m_tuple_size>{});
}

}  // namespace detray