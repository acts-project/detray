/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/// The detray detector volume.
///
/// Volume class that acts as a logical container in the geometry for geometry
/// objects, such as module surfaces or portals. The information is kept as
/// index ranges into larger containers that are owned by the detector
/// implementation. For every object type a different range can be set (e.g
/// for sensitive surfaces or portals).
///
/// @tparam object_registry_t the type of objects contained in the volume
/// @tparam scalar_t type of scalar used in the volume
/// @tparam sf_finder_link_t the type of link to the volumes surfaces finder
///         (geometry surface finder structure). The surface finder types
///         cannot be given directly, since the containers differ between host
///         and device. The surface finders reside in an 'unrollable
///         container' and are called per volume in the navigator during local
///         navigation.
/// @tparam array_t the type of the internal array, must have STL semantics
template <typename object_registry_t, typename scalar_t,
          typename sf_finder_link_t = dindex,
          template <typename, std::size_t> class array_t = darray>
class volume {

    public:
    /// The types of surfaces a volume can contain (modules, portals)
    using objects = object_registry_t;
    /// ID type of a surface
    using obj_link_type = typename objects::link_type;
    /// How to address objects in a container. Can be a range or might become
    /// an index if surfaces are addressed as batches
    using obj_range_type = typename objects::range_type;
    /// Type and index of the surface finder strucure used by this volume
    using sf_finder_link_type = sf_finder_link_t;

    /// In case the geometry needs to be printed
    using name_map = std::map<dindex, std::string>;
    /// Voume tag, used for sfinae
    using volume_def = volume<objects, scalar_t, sf_finder_link_type, array_t>;

    /// Default constructor
    volume() = default;

    /// Constructor from boundary values.
    ///
    /// @param bounds values of volume boundaries. They depend on the volume
    ///               shape, which is defined by its portals and are chosen in
    ///               the detector builder
    volume(const array_t<scalar_t, 6> &bounds) : _bounds(bounds) {}

    /// @return the bounds - const access
    DETRAY_HOST_DEVICE
    auto bounds() const -> const array_t<scalar_t, 6> & { return _bounds; }

    /// @return the name (add an offset for the detector name)
    DETRAY_HOST_DEVICE
    auto name(const name_map &names) const -> const std::string & {
        return names.at(_index + 1);
    }

    /// @return the index of the volume in the detector container
    DETRAY_HOST_DEVICE
    auto index() const -> dindex { return _index; }

    /// @param index the index of the volume in the detector container
    DETRAY_HOST
    auto set_index(const dindex index) -> void { _index = index; }

    /// set surface finder during detector building
    DETRAY_HOST
    auto set_sf_finder(const sf_finder_link_type &link) -> void {
        _sf_finder = link;
    }

    /// set surface finder during detector building
    DETRAY_HOST
    auto set_sf_finder(const typename sf_finder_link_t::id_type id,
                       const typename sf_finder_link_t::index_type index)
        -> void {
        _sf_finder = {id, index};
    }

    /// @return the surface finder link associated with the volume
    DETRAY_HOST_DEVICE
    auto sf_finder_link() const -> const sf_finder_link_type & {
        return _sf_finder;
    }

    /// @return the type of surface finder associated with the volume
    DETRAY_HOST_DEVICE
    auto sf_finder_type() const -> typename sf_finder_link_t::id_type {
        return detail::get<0>(_sf_finder);
    }

    /// @return the index of the surface finder associated with volume
    DETRAY_HOST_DEVICE
    auto sf_finder_index() const -> typename sf_finder_link_t::index_type {
        return detail::get<1>(_sf_finder);
    }

    /// @return if the volume is empty or not
    DETRAY_HOST_DEVICE
    auto empty() const -> bool { return n_objects<objects::e_surface>() == 0; }

    /// @return the number of surfaces in the volume
    template <typename objects::id range_id = objects::e_surface>
    DETRAY_HOST_DEVICE auto n_objects() const -> dindex {
        return n_in_range(range<range_id>());
    }

    /// Set or update the index into a geometry container identified by the
    /// range_id.
    ///
    /// @param other Surface index range
    template <typename objects::id range_id = objects::e_surface>
    DETRAY_HOST auto update_range(const obj_range_type &other) -> void {
        auto &rg = detail::get<range_id>(_ranges);
        // Range not set yet - initialize
        constexpr obj_range_type empty{};
        if (rg == empty) {
            rg = other;
        } else {
            // Update
            assert(detail::get<1>(rg) == detail::get<0>(other));
            detail::get<1>(rg) = detail::get<1>(other);
        }
    }

    /// @return range of surfaces by surface type - const access.
    template <typename object_t>
    DETRAY_HOST_DEVICE inline constexpr auto range() const
        -> const obj_range_type & {
        constexpr auto index = objects::template get_index<object_t>::value;
        return detail::get<index>(_ranges);
    }

    /// @return range of surfaces- const access.
    template <typename objects::id range_id = objects::e_surface>
    DETRAY_HOST_DEVICE inline constexpr auto range() const
        -> const obj_range_type & {
        return detail::get<range_id>(_ranges);
    }

    /// @return all ranges in the volume.
    DETRAY_HOST_DEVICE
    inline constexpr auto ranges() const
        -> const array_t<obj_range_type, objects::n_types> & {
        return _ranges;
    }

    /// @return the number of elements in a given range
    template <typename range_type>
    DETRAY_HOST_DEVICE inline auto n_in_range(range_type &&rg) const -> dindex {
        return detail::get<1>(rg) - detail::get<0>(rg);
    }

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    auto operator==(const volume &rhs) const -> bool {
        return (_bounds == rhs._bounds && _index == rhs._index &&
                _ranges == rhs._ranges && _sf_finder == rhs._sf_finder);
    }

    private:
    /// Bounds section, default for r, z, phi
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
    array_t<obj_range_type, objects::n_types> _ranges = {};

    /// Links to a specific sf_finder structure for this volume
    sf_finder_link_type _sf_finder = {};
};

}  // namespace detray
