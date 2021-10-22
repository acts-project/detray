/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "utils/indexing.hpp"

namespace detray {

/** Volume class that acts as a logical container in the geometry for geometry
 *  objects, such as module surfaces or portals. The information is kept as
 *  index ranges in larger containers that are owned by either the geometry or
 *  the detector implementations.
 *
 * @tparam array_type the type of the internal array, must have STL semantics
 */
template <typename object_ids, typename range_type = dindex_range,
          template <typename, unsigned int> class array_type = darray>
class volume {

    public:
    // The type of objects a volume can contain
    using objects = object_ids;
    // In case the geometry needs to be printed
    using name_map = std::map<dindex, std::string>;

    /** Contructor with bounds
     * @param bounds of the volume
     */
    volume(const array_type<scalar, 6> &bounds) : _bounds(bounds) {}

    /** @return the bounds - const access */
    inline const array_type<scalar, 6> &bounds() const { return _bounds; }

    /** @return the name */
    inline const std::string &name(const name_map &names) const {
        return names.at(_index);
    }

    /** @return the index */
    inline dindex index() const { return _index; }

    /** @param index the index */
    inline void set_index(const dindex index) { _index = index; }

    /** @return the entry into the local surface finders */
    inline dindex surfaces_finder_entry() const {
        return _surfaces_finder_entry;
    }

    /** @param entry the entry into the local surface finders */
    inline void set_surfaces_finder(const dindex entry) {
        _surfaces_finder_entry = entry;
    }

    /** @return if the volume is empty or not */
    inline bool empty() const {
        return n_objects<object_ids::e_surface>() == 0;
    }

    /** @return the number of surfaces in the volume */
    template <object_ids range_id = object_ids::e_surface>
    inline auto n_objects() const {
        return n_in_range(range<range_id>());
    }

    /** Set or update the index into a geometry container identified by the
     *  range_id.
     *
     * @param other Surface index range
     */
    template <object_ids range_id = object_ids::e_surface>
    inline void set_range(const range_type &other) {
        auto &rg = std::get<range_id>(_ranges);
        // Range not set yet - initialize
        constexpr range_type empty{};
        if (rg == empty) {
            rg = other;
            // Update
        } else {
            std::get<1>(rg) += n_in_range(other);
        }
    }

    /** @return range of surfaces- const access */
    template <object_ids range_id = object_ids::e_surface>
    inline constexpr const auto &range() const {
        return std::get<range_id>(_ranges);
    }

    /** @return range of surfaces and portals (must be contiguous!) */
    inline const auto full_range() const {
        // There may be volumes without surfaces, but never without portals
        if ((_surface_range[0] + _surface_range[1] == dindex_invalid) or
            (n_in_range(_surface_range) == 0)) {
            return _portal_range;
        } else if (_portal_range[0] < _surface_range[0]) {
            return dindex_range{_portal_range[0], _surface_range[1]};
        } else {
            return dindex_range{_surface_range[0], _portal_range[1]};
        }
    }

    private:

    /** @return the number of elements in a given range */
    template <typename range_t>
    inline dindex n_in_range(range_t &&rg) const {
        return std::get<1>(rg) - std::get<0>(rg);
    }

    private:
    /** Bounds section, default for r, z, phi */
    array_type<scalar, 6> _bounds = {0.,
                                     std::numeric_limits<scalar>::max(),
                                     -std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::max(),
                                     -M_PI,
                                     M_PI};

    /** Volume index */
    dindex _index = dindex_invalid;

    /** Ranges in geometry containers for different objects types that belong
     * to this volume */
    array_type<dindex_range, object_ids::e_object_types> _ranges = {};

    /** Index into the surface finder container */
    dindex _surfaces_finder_entry = dindex_invalid;
};

}  // namespace detray