/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <map>
#include <string>

#include "utils/indexing.hpp"

namespace detray {

/** Volume class that holds the local information of the volume its surfaces
 * and portals. Everything is kept as index ranges in larger containers that
 * are owned by either the geometry or the detector implementations.
 *
 * @tparam array_type the type of the internal array, must have STL semantics
 */
template <template <typename, unsigned int> class array_type = darray>
class volume {

    public:
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
    inline bool empty() const { return is_empty_range(_surface_range); }

    /** @return the number of surfaces in the volume */
    template <bool primitive = true>
    inline dindex n_objects() {
        if constexpr (primitive) {
            return n_in_range(_surface_range);
        } else {
            return n_in_range(_portal_range);
        }
    }

    /** @return the number of surfaces in the volume */
    template <bool primitive = true>
    inline const dindex n_objects() const {
        if constexpr (primitive) {
            return n_in_range(_surface_range);
        } else {
            return n_in_range(_portal_range);
        }
    }

    /** Set the index into the geometry surface container
     *
     * @param range Surface index range
     */
    template <bool surface_range = true>
    inline void set_range(dindex_range range) {
        if constexpr (surface_range) {
            update_range(_surface_range, std::move(range));
        } else {
            update_range(_portal_range, std::move(range));
        }
    }

    /** @return range of surfaces- const access */
    template <bool surface_range = true>
    inline const auto &range() const {
        if constexpr (surface_range) {
            return _surface_range;
        } else {
            return _portal_range;
        }
    }

    /** Print volume.
     *
     * @returns the volume description as a string
     */
    inline const std::string to_string(const name_map &names) const {
        std::stringstream ss;
        ss << " - name: '" << name(names) << "'" << std::endl;

        ss << "     contains    " << n_objects<true>() << " surfaces "
           << std::endl;

        ss << "                 " << n_objects<false>() << " portals "
           << std::endl;

        if (surfaces_finder_entry() != dindex_invalid) {
            ss << "  sf finders idx " << surfaces_finder_entry() << std::endl;
        }

        ss << "     bounds r = (" << _bounds[0] << ", " << _bounds[1] << ")"
           << std::endl;
        ss << "            z = (" << _bounds[2] << ", " << _bounds[3] << ")"
           << std::endl;

        return ss.str();
    };

    private:
    /**
     * @param range Any index range
     *
     * @return the number of indexed objects
     */
    inline dindex n_in_range(const dindex_range &range) {
        return range[1] - range[0];
    }

    /**
     * @param range Any index range
     *
     * @return the number of indexed objects
     */
    inline const dindex n_in_range(const dindex_range &range) const {
        return range[1] - range[0];
    }

    /** Test whether a range is empty
     *
     * @param range Any index range
     *
     * @return boolean whether the range is empty
     */
    inline bool is_empty_range(const dindex_range &range) {
        return n_in_range(range) == 0;
    }

    /** Test whether a range is empty - const
     *
     * @param range Any index range
     *
     * @return boolean whether the range is empty
     */
    inline const bool is_empty_range(const dindex_range &range) const {
        return n_in_range(range) == 0;
    }

    /** Set or update a range
     *
     * @param range One of the volume member ranges
     * @param other new index range
     *
     * @return boolean whether the range is empty
     */
    inline void update_range(dindex_range &range, dindex_range &&other) {
        // Range not set yet
        if (range[0] == dindex_invalid) {
            range = other;
        } else {
            range[1] += other[1] - other[0];
        }
    }

    /** Bounds section, default for r, z, phi */
    array_type<scalar, 6> _bounds = {0.,
                                     std::numeric_limits<scalar>::max(),
                                     -std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::max(),
                                     -M_PI,
                                     M_PI};

    /** Volume index */
    dindex _index = dindex_invalid;

    /** Index ranges in the detector surface/portal containers.*/
    dindex_range _surface_range = {dindex_invalid, dindex_invalid};
    dindex_range _portal_range = {dindex_invalid, dindex_invalid};

    /** Index into the surface finder container */
    dindex _surfaces_finder_entry = dindex_invalid;
};

}  // namespace detray