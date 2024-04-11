/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/// @brief Wraps a jagged vector-like container and implements a data store
/// around it.
///
/// @tparam T The type of the collection data, e.g. transforms
/// @tparam container_t The type of container to use for the data collection.
/// @tparam context_t the context with which to retrieve the correct data.
template <typename T, template <typename...> class container_t,
          typename context_t>
class context_aware_store {
    public:
    using base_type = container_t<T>;
    using vector_type = typename base_type::value_type;
    using context_type = context_t;

    // TODO: do we need this?
    /*
    /// How to find data in the store
    /// @{
    using link_type = dindex;
    using single_link = dindex;
    using range_link = dindex_range;
    /// @}
    */

    /// Vecmem view types
    using view_type = detail::get_view_t<container_t<T>>;
    /*
     TODO: uncommenting the next line results in compilation problems:
     Originated from container_views.hpp
     error: more than one partial specialization matches the template argument
     list of class "detray::detail::has_view<.......>"
     "detray::detail::has_view<const vecmem::vector<T>, void>"
     "detray::detail::has_view<const vecmem::jagged_vector<T>, void>"
    */
    // using const_view_type = detail::get_view_t<const container_t<T>>;

    /// Empty container
    constexpr context_aware_store() = default;

    // Delegate constructors to container, which handles the memory

    /// Copy construct from element types
    // TODO: do we need this?
    /* constexpr explicit single_store(const T &arg) : m_container(arg) {} */

    /// Construct with a specific memory resource @param resource
    /// (host-side only)
    template <typename allocator_t = vecmem::memory_resource,
              std::enable_if_t<not detail::is_device_view_v<allocator_t>,
                               bool> = true>
    DETRAY_HOST explicit context_aware_store(allocator_t& resource)
        : m_container(&resource) {}

    /// Copy Construct with a specific memory resource @param resource
    /// (host-side only)
    // TODO: do we need this?
    /*
    template <typename allocator_t = vecmem::memory_resource,
              typename C = container_t<T>,
              std::enable_if_t<std::is_same_v<C, std::vector<T>>, bool> = true>
    DETRAY_HOST explicit context_aware_store(allocator_t &resource, const T
    &arg) : m_container(&resource, arg) {}
    */

    /// Construct from the container @param view . Mainly used device-side.
    template <typename container_view_t,
              std::enable_if_t<detail::is_device_view_v<container_view_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE context_aware_store(container_view_t& view)
        : m_container(view) {}

    /// @returns a pointer to the underlying container - const
    DETRAY_HOST_DEVICE
    constexpr auto data() const noexcept -> const base_type* {
        return &m_container;
    }

    /// @returns a pointer to the underlying container - non-const
    DETRAY_HOST_DEVICE
    constexpr auto data() noexcept -> base_type* { return &m_container; }

    /// @returns the size of the underlying container
    DETRAY_HOST_DEVICE
    constexpr auto size(const context_type& ctx) const noexcept -> dindex {
        return static_cast<dindex>(m_container.at(ctx.get()).size());
    }

    /// @returns true if the underlying container is empty
    DETRAY_HOST_DEVICE
    constexpr auto empty(const context_type& ctx) const noexcept -> bool {
        return m_container.at(ctx.get()).empty();
    }

    /// @returns the collections iterator at the start position.
    DETRAY_HOST_DEVICE
    constexpr auto begin(const context_type& ctx) {
        return m_container.at(ctx.get()).begin();
    }

    /// @returns the collections iterator sentinel.
    DETRAY_HOST_DEVICE
    constexpr auto end(const context_type& ctx) {
        return m_container.at(ctx.get()).end();
    }

    /// @returns access to the second-level vector in underlying container -
    /// const
    DETRAY_HOST_DEVICE
    auto get(const context_type& ctx) const noexcept -> const vector_type& {
        return m_container.at(ctx.get());
    }

    /// @returns access to the second-level vector in underlying container -
    /// non-const
    DETRAY_HOST_DEVICE
    auto get(const context_type& ctx) noexcept -> vector_type& {
        return m_container.at(ctx.get());
    }

    /// @returns access to the underlying container - const
    DETRAY_HOST_DEVICE
    auto get() const noexcept -> const base_type& { return m_container; }

    /// @returns access to the underlying container - non-const
    DETRAY_HOST_DEVICE
    auto get() noexcept -> base_type& { return m_container; }

    // TODO: does it make sense to introduce the operator[] for this store?
    /// Elementwise access. Needs @c operator[] for storage type - non-const
    /*
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const dindex i) {
      return m_container[i];
    }

    /// Elementwise access. Needs @c operator[] for storage type - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const dindex i) const {
      return m_container[i];
    }
    */

    /// @returns context based access to an element (also range checked)
    DETRAY_HOST_DEVICE
    constexpr auto at(const dindex i, const context_type& ctx) const noexcept
        -> const T& {
        return m_container.at(ctx.get()).at(i);
    }

    /// @returns context based access to an element (also range checked)
    DETRAY_HOST_DEVICE
    constexpr auto at(const dindex i, const context_type& ctx) noexcept -> T& {
        return m_container.at(ctx.get()).at(i);
    }

    /// Removes and destructs all elements in the container for a given context
    DETRAY_HOST void clear(const context_type& ctx) {
        m_container.at(ctx.get()).clear();
    }

    /// Reserve memory of size @param n for a given geometry context
    DETRAY_HOST void reserve(std::size_t n, const context_type& ctx) {
        if (ctx.get() <= m_container.size()) {
            auto& vect =
                ctx.get() == m_container.size()
                    ? *(m_container.insert(m_container.end(), vector_type()))
                    : m_container.at(ctx.get());
            vect.reserve(n);
        } else {
            // TODO: how to react to these situations?
        }
    }

    /// Resize the underlying container to @param n for a given geometry context
    DETRAY_HOST void resize(std::size_t n, const context_type& ctx) {
        if (ctx.get() <= m_container.size()) {
            auto& vect =
                ctx.get() == m_container.size()
                    ? *(m_container.insert(m_container.end(), vector_type()))
                    : m_container.at(ctx.get());
            vect.resize(n);
        } else {
            // TODO: how to react to these situations?
        }
    }

    /// Add a new element to the collection - copy
    ///
    /// @tparam U type that can be converted to T
    ///
    /// @param arg the constructor argument
    ///
    /// @note in general can throw an exception
    template <typename U>
    DETRAY_HOST constexpr auto push_back(
        const U& arg, const context_type& ctx) noexcept(false) -> void {
        if (ctx.get() <= m_container.size()) {
            auto& vect =
                ctx.get() == m_container.size()
                    ? *(m_container.insert(m_container.end(), vector_type()))
                    : m_container.at(ctx.get());
            vect.push_back(arg);
        } else {
            // TODO: how to react to these situations?
        }
    }

    /// Add a new element to the collection - move
    ///
    /// @tparam U type that can be converted to T
    ///
    /// @param arg the constructor argument
    ///
    /// @note in general can throw an exception
    template <typename U>
    DETRAY_HOST constexpr auto push_back(
        U&& arg, const context_type& ctx) noexcept(false) -> void {
        if (ctx.get() <= m_container.size()) {
            auto& vect =
                ctx.get() == m_container.size()
                    ? *(m_container.insert(m_container.end(), vector_type()))
                    : m_container.at(ctx.get());
            vect.push_back(std::move(arg));
        } else {
            // TODO: how to react to these situations?
        }
    }

    /// Add a new element to the collection in place
    ///
    /// @tparam Args are the types of the constructor arguments
    ///
    /// @param args is the list of constructor arguments
    ///
    /// @note in general can throw an exception
    template <typename... Args>
    DETRAY_HOST constexpr decltype(auto) emplace_back(
        const context_type& ctx, Args&&... args) noexcept(false) {
        if (ctx.get() <= m_container.size()) {
            auto& vect =
                ctx.get() == m_container.size()
                    ? *(m_container.insert(m_container.end(), vector_type()))
                    : m_container.at(ctx.get());
            return vect.emplace_back(std::forward<Args>(args)...);
        } else {
            // TODO: how to react to these situations? For now throw an
            // exception to avoid warnings
            throw std::runtime_error{""};
        }
    }

    /// Insert another collection - copy
    ///
    /// @tparam U type that represents second level vector in a jagged vector
    ///
    /// @param new_data is the new collection to be added for the given context
    ///
    /// @note in general can throw an exception
    template <typename U>
    DETRAY_HOST auto insert(U& new_data,
                            const context_type& ctx) noexcept(false) -> void {
        if (ctx.get() <= m_container.size()) {
            auto& vect =
                ctx.get() == m_container.size()
                    ? *(m_container.insert(m_container.end(), vector_type()))
                    : m_container.at(ctx.get());
            vect.reserve(vect.size() + new_data.size());
            vect.insert(vect.end(), new_data.begin(), new_data.end());
        } else {
            // TODO: how to react to these situations?
        }
    }

    /// Insert another collection - move
    ///
    /// @tparam U type that represents second level vector in a jagged vector
    ///
    /// @param new_data is the new collection to be added for the given context
    ///
    /// @note in general can throw an exception
    template <typename U>
    DETRAY_HOST auto insert(U&& new_data,
                            const context_type& ctx) noexcept(false) -> void {
        if (ctx.get() <= m_container.size()) {
            auto& vect =
                ctx.get() == m_container.size()
                    ? *(m_container.insert(m_container.end(), vector_type()))
                    : m_container.at(ctx.get());
            vect.reserve(vect.size() + new_data.size());
            vect.insert(vect.end(), std::make_move_iterator(new_data.begin()),
                        std::make_move_iterator(new_data.end()));
        } else {
            // TODO: how to react to these situations?
        }
    }

    /// Append the content from another store to the current one for a given
    /// context
    ///
    /// @param other The other container
    ///
    /// @note in general can throw an exception
    DETRAY_HOST void append(context_aware_store& other,
                            const context_type& ctx) noexcept(false) {
        insert(other.m_container.at(ctx.get()), ctx);
    }

    /// Append the content from another store to the current one for a given
    /// context - move
    ///
    /// @param other The other container
    ///
    /// @note in general can throw an exception
    DETRAY_HOST void append(context_aware_store&& other,
                            const context_type& ctx) noexcept(false) {
        insert(std::move(other.m_container.at(ctx.get())), ctx);
    }

    /// @return the view on the underlying container - non-const
    DETRAY_HOST auto get_data() -> view_type {
        return detray::get_data(m_container);
    }

    // TODO: uncomment this after the issue with const_view_type gets sorted out
    /*
    /// @return the view on the underlying container - const
    DETRAY_HOST auto get_data() const -> const_view_type {
        return detray::get_data(m_container);
    }
    */

    private:
    base_type m_container;
};

}  // namespace detray
