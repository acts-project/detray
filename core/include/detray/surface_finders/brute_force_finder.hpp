/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s).
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/// @brief A collection of brute force surface finders, callable by index.
///
/// This class fulfills all criteria to be used in the detector @c multi_store .
///
/// @tparam surface_t the type of surface data handles in a detector.
/// @tparam container_t the types of underlying containers to be used.
template <class surface_t, typename container_t = host_container_types>
class brute_force_collection {

    public:
    template <typename T>
    using vector_type = typename container_t::template vector_type<T>;
    using size_type = dindex;

    /// A nested surface finder that returns all surfaces in a range (brute
    /// force). This type will be returned when the surface collection is
    /// queried for the surfaces of a particular volume.
    struct brute_forcer
        : public detray::ranges::subrange<const vector_type<surface_t>> {

        using base = detray::ranges::subrange<const vector_type<surface_t>>;

        /// Default constructor
        brute_forcer() = default;

        /// Constructor from @param surface_range - move
        DETRAY_HOST_DEVICE constexpr brute_forcer(
            const vector_type<surface_t>& surfaces, const dindex_range& range)
            : base(surfaces, range) {}

        /// @returns the complete surface range of the search volume
        template <typename detector_t, typename track_t>
        DETRAY_HOST_DEVICE constexpr auto search(
            const detector_t& /*det*/,
            const typename detector_t::volume_type& /*volume*/,
            const track_t& /*track*/) const {
            return *this;
        }

        /// @returns an iterator over all surfaces in the data structure
        DETRAY_HOST_DEVICE constexpr auto all() const { return *this; }

        /// @return the maximum number of surface candidates during a
        /// neighborhood lookup
        DETRAY_HOST_DEVICE constexpr auto n_max_candidates() const
            -> unsigned int {
            return static_cast<unsigned int>(this->size());
        }
    };

    using value_type = brute_forcer;

    using view_type =
        dmulti_view<dvector_view<size_type>, dvector_view<surface_t>>;
    using const_view_type = dmulti_view<dvector_view<const size_type>,
                                        dvector_view<const surface_t>>;

    /// Default constructor
    constexpr brute_force_collection() {
        // Start of first subrange
        m_offsets.push_back(0u);
    };

    /// Constructor from memory resource
    DETRAY_HOST
    explicit constexpr brute_force_collection(vecmem::memory_resource* resource)
        : m_offsets(resource), m_surfaces(resource) {
        // Start of first subrange
        m_offsets.push_back(0u);
    }

    /// Device-side construction from a vecmem based view type
    template <typename coll_view_t,
              typename std::enable_if_t<detail::is_device_view_v<coll_view_t>,
                                        bool> = true>
    DETRAY_HOST_DEVICE brute_force_collection(coll_view_t& view)
        : m_offsets(detail::get<0>(view.m_view)),
          m_surfaces(detail::get<1>(view.m_view)) {}

    /// @returns number of surface collections (at least on per volume) - const
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> size_type {
        // The start index of the first range is always present
        return static_cast<dindex>(m_offsets.size()) - 1u;
    }

    /// @note outside of navigation, the number of elements is unknown
    DETRAY_HOST_DEVICE
    constexpr auto empty() const noexcept -> bool {
        return size() == size_type{0};
    }

    /// @return access to the surface container - const.
    DETRAY_HOST_DEVICE
    auto all() const -> const vector_type<surface_t>& { return m_surfaces; }

    /// @return access to the surface container - non-const.
    DETRAY_HOST_DEVICE
    auto all() -> vector_type<surface_t>& { return m_surfaces; }

    /// Create brute force surface finder from surface container - const
    DETRAY_HOST_DEVICE
    auto operator[](const size_type i) const -> value_type {
        return {m_surfaces, dindex_range{m_offsets[i], m_offsets[i + 1u]}};
    }

    /// Add a new surface collection
    template <
        typename sf_container_t,
        typename std::enable_if_t<detray::ranges::range_v<sf_container_t>,
                                  bool> = true,
        typename std::enable_if_t<
            std::is_same_v<typename sf_container_t::value_type, surface_t>,
            bool> = true>
    DETRAY_HOST auto push_back(const sf_container_t& surfaces) noexcept(false)
        -> void {
        m_surfaces.reserve(m_surfaces.size() + surfaces.size());
        m_surfaces.insert(m_surfaces.end(), surfaces.begin(), surfaces.end());
        // End of this range is the start of the next range
        m_offsets.push_back(static_cast<dindex>(m_surfaces.size()));
    }

    /// @return the view on the brute force finders - non-const
    DETRAY_HOST
    constexpr auto get_data() noexcept -> view_type {
        return {detray::get_data(m_offsets), detray::get_data(m_surfaces)};
    }

    /// @return the view on the brute force finders - const
    DETRAY_HOST
    constexpr auto get_data() const noexcept -> const_view_type {
        return {detray::get_data(m_offsets), detray::get_data(m_surfaces)};
    }

    private:
    /// Offsets for the respective volumes into the surface storage
    vector_type<size_type> m_offsets{};
    /// The storage for all surface handles
    vector_type<surface_t> m_surfaces{};
};

}  // namespace detray
