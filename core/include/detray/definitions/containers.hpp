/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple.hpp"

// Vecmem include(s)
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/containers/vector.hpp>

// System include(s)
#include <map>
#include <type_traits>
#include <vector>

namespace detray {

template <typename value_t, std::size_t kDIM>
using darray = std::array<value_t, kDIM>;

template <typename value_t>
using dvector = vecmem::vector<value_t>;

template <typename value_t>
using djagged_vector = vecmem::jagged_vector<value_t>;

template <typename key_t, typename value_t>
using dmap = std::map<key_t, value_t>;

template <class... types>
using dtuple = detray::tuple<types...>;

template <typename value_t>
struct compact_device_vector {
    using value_type = value_t;
    using size_type = unsigned int;
    using difference_type = std::ptrdiff_t;
    using reference = std::add_lvalue_reference_t<value_type>;
    using const_reference =
        std::add_lvalue_reference_t<std::add_const_t<value_type>>;
    using pointer = std::add_pointer_t<value_type>;
    using const_pointer = std::add_pointer_t<std::add_const_t<value_type>>;

    using iterator = pointer;
    using const_iterator = const_pointer;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    DETRAY_HOST_DEVICE compact_device_vector(
        const vecmem::data::vector_view<value_type>& data)
        : m_ptr(data.ptr()), m_size(data.size()) {}

    DETRAY_HOST_DEVICE
    reference at(size_type pos) {
        assert(pos < size());
        return m_ptr[pos];
    }
    DETRAY_HOST_DEVICE
    const_reference at(size_type pos) const {
        assert(pos < size());
        return m_ptr[pos];
    }

    DETRAY_HOST_DEVICE
    reference operator[](size_type pos) { return m_ptr[pos]; }
    DETRAY_HOST_DEVICE
    const_reference operator[](size_type pos) const { return m_ptr[pos]; }

    DETRAY_HOST_DEVICE
    reference front() {
        assert(size() > 0);
        return m_ptr[0];
    }
    DETRAY_HOST_DEVICE
    const_reference front() const {
        assert(size() > 0);
        return m_ptr[0];
    }

    DETRAY_HOST_DEVICE
    reference back() {
        assert(size() > 0);
        return m_ptr[size() - 1];
    }
    DETRAY_HOST_DEVICE
    const_reference back() const {
        assert(size() > 0);
        return m_ptr[size() - 1];
    }

    DETRAY_HOST_DEVICE
    pointer data() { return m_ptr; }
    DETRAY_HOST_DEVICE
    const_pointer data() const { return m_ptr; }

    DETRAY_HOST_DEVICE
    iterator begin() { return iterator(m_ptr); }
    DETRAY_HOST_DEVICE
    const_iterator begin() const { return const_iterator(m_ptr); }
    DETRAY_HOST_DEVICE
    const_iterator cbegin() const { return begin(); }

    DETRAY_HOST_DEVICE
    iterator end() { return iterator(m_ptr + size()); }
    DETRAY_HOST_DEVICE
    const_iterator end() const { return const_iterator(m_ptr + size()); }
    DETRAY_HOST_DEVICE
    const_iterator cend() const { return end(); }

    DETRAY_HOST_DEVICE
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    DETRAY_HOST_DEVICE
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end());
    }
    DETRAY_HOST_DEVICE
    const_reverse_iterator crbegin() const { return rbegin(); }

    DETRAY_HOST_DEVICE
    reverse_iterator rend() { return reverse_iterator(begin()); }
    DETRAY_HOST_DEVICE
    const_reverse_iterator rend() const {
        return const_reverse_iterator(begin());
    }
    DETRAY_HOST_DEVICE
    const_reverse_iterator crend() const { return rend(); }

    DETRAY_HOST_DEVICE
    bool empty() const { return (size() == 0); }
    DETRAY_HOST_DEVICE
    size_type size() const { return m_size; }
    DETRAY_HOST_DEVICE
    size_type capacity() const { return m_size; }

    private:
    pointer m_ptr;
    size_type m_size;
};

/// @brief Bundle container type definitions
template <template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector,
          template <typename, typename> class map_t = dmap>
struct container_types {
    template <typename T>
    using vector_type = vector_t<T>;

    template <class... T>
    using tuple_type = dtuple<T...>;

    template <typename T, std::size_t kDIM>
    using array_type = darray<T, kDIM>;

    template <typename T>
    using jagged_vector_type = jagged_vector_t<T>;

    template <typename K, typename T>
    using map_type = map_t<K, T>;
};

/// Defining some common types
using host_container_types = container_types<>;

namespace detail {

// make std::get available in detray detail namespace, where also the thrust and
// index specific overloads live.
using std::get;

/// Trait class to figure out if a given type has a @c reserve(...) function
template <typename T>
struct has_reserve {

    private:
    /// Function returning @c std::true_type for types that do have a @c
    /// reserve(...) function
    template <typename C>
    static constexpr auto check(C*) ->
        typename std::is_void<decltype(std::declval<C>().reserve(
            std::declval<typename C::size_type>()))>::type;

    /// Function returning @c std::false_type for types that fair the previous
    /// function
    template <typename>
    static constexpr std::false_type check(...);

    /// Declare the value type of this trait class
    using type = decltype(check<T>(nullptr));

    public:
    /// Value of the check
    static constexpr bool value = type::value;
};

/// @name Functions calling or not calling reserve(...) based on whether it's
/// available
/// @{
template <typename T>
    requires has_reserve<T>::value
DETRAY_HOST_DEVICE void call_reserve(T& obj, std::size_t newsize) {
    obj.reserve(newsize);
}

template <typename T>
    requires(!has_reserve<T>::value)
DETRAY_HOST_DEVICE void call_reserve(T& /*obj*/, std::size_t /*newsize*/) {
    /*Not defined*/
}
/// @}

}  // namespace detail

}  // namespace detray
