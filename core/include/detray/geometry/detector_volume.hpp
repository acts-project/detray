/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// @brief The detray detector volume.
///
/// Volume class that acts as a logical container in the detector for geometry
/// objects, i.e. surfaces. The volume boundary surfaces, the so-called portals,
/// carry index links that join adjacent volumes. The volume class itself does
/// not contain any data itself, but keeps index-based links into the data
/// containers that are managed by the detector.
/// Every type of surface that is known by the volume (determined by the type
/// and size of the ID enum) lives in it's own geometry accelerator data
/// structure, e.g. portals reside in a brute force accelerator (a simple
/// vector), while sensitive surfaces are usually sorted into a spacial grid.
///
/// @tparam ID enum of object types contained in the volume
///         (@see @c detector_metadata ).
/// @tparam link_t the type of link to the volumes surfaces finder(s)
///         (accelerator structure, e.g. a grid). The surface finder types
///         cannot be given directly, since the containers differ between host
///         and device. The surface finders reside in an 'unrollable
///         container' and are called per volume in the navigator during local
///         navigation.
/// @tparam scalar_t type of scalar used in the volume.
/// @tparam array_t the type of the internal array, must have STL semantics.
template <typename ID, typename link_t = dtyped_index<dindex, dindex>,
          typename scalar_t = scalar,
          template <typename, std::size_t> class array_t = darray>
class detector_volume {

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
    /// portals (the accelerator structures id and index).
    using link_type = dmulti_index<link_t, ID::e_size>;

    /// In case the geometry needs to be printed
    using name_map = std::map<dindex, std::string>;

    /// Voume tag, used for sfinae
    using volume_def = detector_volume<ID, link_t, scalar_t, array_t>;

    /// Default constructor builds an infinitely long cylinder
    constexpr detector_volume() = default;

    /// Constructor from boundary values.
    ///
    /// @param id id values that determines how to interpret the bounds.
    /// @param bounds values of volume boundaries. They depend on the volume
    ///               shape, which is defined by its portals and are chosen in
    ///               the detector builder.
    constexpr detector_volume(const volume_id id,
                              const array_t<scalar_t, 6> &bounds)
        : _id{id}, _bounds(bounds) {}

    /// @return the volume shape id, e.g. 'cylinder'
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> volume_id { return _id; }

    /// @return the bounds - const access
    DETRAY_HOST_DEVICE
    constexpr auto bounds() const -> const array_t<scalar_t, 6> & {
        return _bounds;
    }

    /// @return the volume name (add an offset for the detector name)
    DETRAY_HOST_DEVICE
    constexpr auto name(const name_map &names) const -> const std::string & {
        return names.at(_index + 1);
    }

    /// @return the index of the volume in the detector container
    DETRAY_HOST_DEVICE
    constexpr auto index() const -> dindex { return _index; }

    /// @param index the index of the volume in the detector container
    DETRAY_HOST
    constexpr auto set_index(const dindex index) -> void { _index = index; }

    /// @return link of a type of object - const access.
    template <ID obj_id = ID::e_sensitive>
    DETRAY_HOST_DEVICE constexpr auto link() const -> const link_t & {
        return detail::get<obj_id>(_sf_finder_links);
    }

    /// @return link of a type of object - const access.
    template <ID obj_id = ID::e_sensitive>
    DETRAY_HOST_DEVICE constexpr auto link() -> link_t & {
        return detail::get<obj_id>(_sf_finder_links);
    }

    /// set surface finder during detector building
    template <ID obj_id = ID::e_sensitive>
    DETRAY_HOST constexpr auto set_link(const link_t &link) -> void {
        _sf_finder_links[obj_id] = link;
    }

    /// set surface finder during detector building
    template <ID obj_id = ID::e_sensitive>
    DETRAY_HOST constexpr auto set_link(const typename link_t::id_type id,
                                        const typename link_t::index_type index)
        -> void {
        _sf_finder_links[obj_id] = link_t{id, index};
    }

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const detector_volume &rhs) const -> bool {
        return (_bounds == rhs._bounds && _index == rhs._index &&
                _sf_finder_links == rhs._sf_finder_links &&
                _sf_finder_links[ID::e_sensitive] ==
                    rhs._sf_finder_links[ID::e_sensitive]);
    }

    private:
    /// How to interpret the boundary values
    volume_id _id = volume_id::e_cylinder;

    /// Bounds section, default for cylinder volume
    array_t<scalar_t, 6> _bounds = {0.,
                                    std::numeric_limits<scalar_t>::max(),
                                    -std::numeric_limits<scalar_t>::max(),
                                    std::numeric_limits<scalar_t>::max(),
                                    -M_PI,
                                    M_PI};

    /// Volume index in the detector's volume container
    dindex _index = dindex_invalid;

    /// Indices in geometry containers for different objects types are
    /// contained in this volume
    link_type _sf_finder_links = {};
};

}  // namespace detray
