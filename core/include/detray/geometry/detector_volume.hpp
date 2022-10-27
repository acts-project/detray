/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// The detray detector volume.
///
/// Volume class that acts as a logical container in the detector for geometry
/// objects, such as sensitive surfaces or portals. The information is kept as
/// index links into larger containers that are owned by the detector
/// implementation. For every object type a different link can be set (e.g
/// for sensitive surfaces or portals).
/// Furthermore, volumes are associated with geometry accelerator structure,
/// like e.g. a spacial grid, that is used during navigation.
///
/// @tparam ID enum of type of objects contained in the volume
///         (@see @c detector_metadata ).
/// @tparam obj_link_t the type (and number) of links to geometry objects.
/// @tparam sf_finder_link_t the type of link to the volumes surfaces finder(s)
///         (geometry surface finder structure). The surface finder types
///         cannot be given directly, since the containers differ between host
///         and device. The surface finders reside in an 'unrollable
///         container' and are called per volume in the navigator during local
///         navigation.
/// @tparam scalar_t type of scalar used in the volume.
/// @tparam array_t the type of the internal array, must have STL semantics.
template <typename ID, typename obj_link_t = dmulti_index<dindex_range, 1>,
          typename sf_finder_link_t = dtyped_index<dindex, dindex>,
          typename scalar_t = scalar,
          template <typename, std::size_t> class array_t = darray>
class detector_volume {

    public:
    /// How to address objects in a container. Can be a range or might become
    /// an index if surfaces are addressed as batches
    using obj_link_type = obj_link_t;
    /// Type and index of the surface finder strucure used by this volume
    using sf_finder_link_type = sf_finder_link_t;

    /// In case the geometry needs to be printed
    using name_map = std::map<dindex, std::string>;
    /// Voume tag, used for sfinae
    using volume_def = detector_volume<ID, obj_link_type, sf_finder_link_type,
                                       scalar_t, array_t>;

    /// Default constructor
    constexpr detector_volume() = default;

    /// Constructor from boundary values.
    ///
    /// @param bounds values of volume boundaries. They depend on the volume
    ///               shape, which is defined by its portals and are chosen in
    ///               the detector builder
    constexpr detector_volume(const array_t<scalar_t, 6> &bounds)
        : _bounds(bounds) {}

    /// @return the bounds - const access
    DETRAY_HOST_DEVICE
    constexpr auto bounds() const -> const array_t<scalar_t, 6> & {
        return _bounds;
    }

    /// @return the name (add an offset for the detector name)
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

    /// set surface finder during detector building
    DETRAY_HOST
    constexpr auto set_sf_finder(const sf_finder_link_type &link) -> void {
        _sf_finder = link;
    }

    /// set surface finder during detector building
    DETRAY_HOST
    constexpr auto set_sf_finder(
        const typename sf_finder_link_t::id_type id,
        const typename sf_finder_link_t::index_type index) -> void {
        _sf_finder = {id, index};
    }

    /// @return the surface finder link associated with the volume
    DETRAY_HOST_DEVICE
    constexpr auto sf_finder_link() const -> const sf_finder_link_type & {
        return _sf_finder;
    }

    /// @return the type of surface finder associated with the volume
    DETRAY_HOST_DEVICE
    constexpr auto sf_finder_type() const ->
        typename sf_finder_link_t::id_type {
        return detail::get<0>(_sf_finder);
    }

    /// @return the index of the surface finder associated with volume
    DETRAY_HOST_DEVICE
    constexpr auto sf_finder_index() const ->
        typename sf_finder_link_t::index_type {
        return detail::get<1>(_sf_finder);
    }

    /// @return link of a type of object - const access.
    template <ID obj_id = ID::e_sensitive>
    DETRAY_HOST_DEVICE constexpr auto obj_link() const -> const
        typename obj_link_type::index_type & {
        return detail::get<obj_id>(_obj_links);
    }

    /// @return link of a type of object - const access.
    template <ID obj_id = ID::e_sensitive>
    DETRAY_HOST_DEVICE constexpr auto obj_link() ->
        typename obj_link_type::index_type & {
        return detail::get<obj_id>(_obj_links);
    }

    /// @returns the entire span of all links combined.
    DETRAY_HOST_DEVICE
    constexpr auto full_range() const -> typename obj_link_type::index_type {
        return {
            detail::get<0>(detail::get<0>(_obj_links)),
            detail::get<1>(detail::get<obj_link_type::size() - 1>(_obj_links))};
    }

    /// @return if the volume is empty or not
    DETRAY_HOST_DEVICE
    constexpr auto empty() const -> bool {
        return n_objects<ID::e_sensitive>() == 0;
    }

    /// @return the number of surfaces in the volume
    /// @note the id has to match the number of index linkss in the volume.
    template <ID obj_id = ID::e_all>
    DETRAY_HOST_DEVICE constexpr auto n_objects() const -> std::size_t {
        if constexpr (obj_id == ID::e_all) {
            return detail::get<1>(
                       detail::get<obj_link_type::size() - 1>(_obj_links)) -
                   detail::get<0>(detail::get<0>(_obj_links));
        }
        return detail::get<1>(obj_link<obj_id>()) -
               detail::get<0>(obj_link<obj_id>());
    }

    /// Set or update the index into a geometry container identified by the
    /// obj_id.
    ///
    /// @note There is no check of overlapping index ranges between the object
    /// types. Use with care!
    ///
    /// @param other Surface index range
    template <ID obj_id = ID::e_sensitive,
              std::enable_if_t<obj_id != ID::e_all, bool> = true>
    DETRAY_HOST auto update_obj_link(
        const typename obj_link_type::index_type &other) noexcept -> void {
        auto &rg = obj_link<obj_id>();
        // Range not set yet - initialize
        constexpr typename obj_link_type::index_type empty{};
        if (rg == empty) {
            rg = other;
        } else {
            // Update
            assert(detail::get<1>(rg) == detail::get<0>(other));
            detail::get<1>(rg) = detail::get<1>(other);
        }
    }

    /// @return all links in the volume.
    DETRAY_HOST_DEVICE
    constexpr auto obj_links() const -> const obj_link_type & {
        return _obj_links;
    }

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const detector_volume &rhs) const -> bool {
        return (_bounds == rhs._bounds && _index == rhs._index &&
                _obj_links == rhs._obj_links && _sf_finder == rhs._sf_finder);
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
    obj_link_type _obj_links = {};

    /// Links to a specific sf_finder structure for this volume
    sf_finder_link_type _sf_finder = {};
};

}  // namespace detray
