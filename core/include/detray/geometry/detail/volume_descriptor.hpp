/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"

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
    explicit constexpr volume_descriptor(const volume_id id) : _id{id} {}

    /// @returns the volume shape id, e.g. 'cylinder'
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> volume_id { return _id; }

    /// @return the index of the volume in the detector volume container.
    DETRAY_HOST_DEVICE
    constexpr auto index() const -> dindex { return _index; }

    /// @param index the index of the volume in the detector volume container.
    DETRAY_HOST
    constexpr auto set_index(const dindex index) -> void { _index = index; }

    /// @returns link to all acceleration data structures - const access
    DETRAY_HOST_DEVICE constexpr auto full_link() const -> const link_type & {
        return _sf_finder_links;
    }

    /// @returns link for a specific type of object - const access.
    template <ID obj_id>
    DETRAY_HOST_DEVICE constexpr auto link() const -> const link_t & {
        return detail::get<obj_id>(_sf_finder_links);
    }

    /// @returns link for a specific type of object - non-const access.
    template <ID obj_id>
    DETRAY_HOST_DEVICE constexpr auto link() -> link_t & {
        return detail::get<obj_id>(_sf_finder_links);
    }

    /// Set surface finder during detector building
    template <ID obj_id>
    DETRAY_HOST constexpr auto set_link(const link_t &link) -> void {
        _sf_finder_links[obj_id] = link;
    }

    /// sSt surface finder during detector building
    template <ID obj_id>
    DETRAY_HOST constexpr auto set_link(const typename link_t::id_type id,
                                        const typename link_t::index_type index)
        -> void {
        _sf_finder_links[obj_id] = link_t{id, index};
    }

    /// Equality operator
    ///
    /// @param rhs is the right-hand side to compare against.
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const volume_descriptor &rhs) const -> bool {
        return (_index == rhs._index &&
                _sf_finder_links == rhs._sf_finder_links &&
                _sf_finder_links[ID::e_sensitive] ==
                    rhs._sf_finder_links[ID::e_sensitive]);
    }

    private:
    /// How to interpret the boundary values
    volume_id _id = volume_id::e_cylinder;

    /// Volume index in the detector's volume container
    dindex _index = dindex_invalid;

    /// Links for every object type to an acceleration data structure
    link_type _sf_finder_links = {};
};

}  // namespace detray
