/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/// @brief The detray detector volume descriptor.
///
/// Contains the data and links to describe a detector volume. This is the type
/// that is stored in the detector data stores.
///
/// @tparam ID enum of object types contained in the volume
///         (@see @c detector_metadata ).
/// @tparam link_t the type of link to the volumes surfaces finder(s)
///         (accelerator structure, e.g. a grid). The surface finder types
///         cannot be given directly, since the containers differ between host
///         and device. The surface finders reside in an 'unrollable tuple
///         container' and are called per volume in the navigator during local
///         navigation.
template <typename ID, typename link_t = dtyped_index<dindex, dindex>>
class volume_descriptor {

    public:
    /// Ids of objects that can be distinguished by the volume
    using object_id = ID;

    /// How to access objects (e.g. sensitives/passives/portals) in this
    /// volume. Keeps one accelerator structure link per object type (by ID):
    ///
    /// link_t : id and index of the accelerator structure in the detector's
    ///          surface store.
    ///
    /// E.g. a 'portal' can be found under @c ID::e_portal in this link,
    /// and will then receive link to the @c brute_force_finder that holds the
    /// portals (the accelerator structure's id and index).
    using link_type = dmulti_index<link_t, ID::e_size>;

    /// Default constructor builds an ~infinitely long cylinder
    constexpr volume_descriptor() = default;

    /// Constructor from shape id.
    ///
    /// @param id id values that determines how to interpret the bounds.
    explicit constexpr volume_descriptor(const volume_id id) : m_id{id} {}

    /// @returns the volume shape id, e.g. 'cylinder'
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> volume_id { return m_id; }

    /// @returns the index of the volume in the detector volume container.
    DETRAY_HOST_DEVICE
    constexpr auto index() const -> dindex { return m_index; }

    /// @param index the index of the volume in the detector volume container.
    DETRAY_HOST
    constexpr auto set_index(const dindex index) -> void { m_index = index; }

    /// @returns the index of the volume tranform in the transform store.
    DETRAY_HOST_DEVICE
    constexpr auto transform() const -> dindex { return m_transform; }

    /// @param index the index of the volume in the detector volume container.
    DETRAY_HOST
    constexpr auto set_transform(const dindex trf_idx) -> void {
        m_transform = trf_idx;
    }

    /// @returns link to all acceleration data structures - const access
    DETRAY_HOST_DEVICE constexpr auto full_link() const -> const link_type & {
        return m_sf_finder_links;
    }

    /// @returns acc data structure link for a specific type of object - const
    template <ID obj_id>
    DETRAY_HOST_DEVICE constexpr auto link() const -> const link_t & {
        return detail::get<obj_id>(m_sf_finder_links);
    }

    /// Set surface finder link from @param link
    template <ID obj_id>
    DETRAY_HOST constexpr auto set_link(const link_t &link) -> void {
        m_sf_finder_links[obj_id] = link;
    }

    /// Set surface finder link from @param id and @param index of the
    /// acceleration data structure (e.g. type and index of grid in surface
    /// store)
    template <ID obj_id>
    DETRAY_HOST constexpr auto set_link(const typename link_t::id_type id,
                                        const typename link_t::index_type index)
        -> void {
        m_sf_finder_links[obj_id] = link_t{id, index};
    }

    /// Equality operator
    ///
    /// @param rhs is the right-hand side to compare against.
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const volume_descriptor &rhs) const -> bool {
        return (m_id == rhs.m_id && m_index == rhs.m_index &&
                m_sf_finder_links == rhs.m_sf_finder_links);
    }

    private:
    /// How to interpret the boundary values
    volume_id m_id = volume_id::e_cylinder;

    /// Volume index in the detector's volume container
    dindex m_index = dindex_invalid;

    /// Volume index in the detector's volume container
    dindex m_transform = dindex_invalid;

    /// Links for every object type to an acceleration data structure
    link_type m_sf_finder_links = {};
};

}  // namespace detray
