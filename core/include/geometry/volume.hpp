/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <string>

namespace detray {

/** Volume class that holds the local information of the volume its surfaces
 * and portals. Everything is kept as index ranges in larger containers that are
 * owned by either the geometry or the detector implementations.
 *
 * @tparam array_type the type of the internal array, must have STL semantics
 */
template <template <typename, unsigned int> class array_type = darray>
class volume {

    public:
    /** Deleted constructor */
    volume() = delete;

    /** Allowed constructors
     * @param name of the volume
     * @param d detector the volume belongs to
     *
     * @note will be contructed boundless
     */
    volume(const std::string &name) : _name(name) {}

    /** Contructor with name and bounds
     * @param name of the volume
     * @param bounds of the volume
     * @param d detector the volume belongs to
     */
    volume(const std::string &name, const array_type<scalar, 6> &bounds)
        : _name(name), _bounds(bounds) {}

    /** Copy ctor makes sure constituents keep valid volume pointer
     *
     * @param other Volume to be copied
     */
    volume(const volume &other) = default;

    /** Equality operator of volumes, convenience function - const
     *
     * @param rhs is the volume to be compared with
     *
     * checks identity for all index ranges and bounds
     **/
    const bool operator==(const volume &rhs) const {
        const bool is_eq_bounds = (_bounds == rhs._bounds);
        const bool is_eq_index = (_index == rhs._index);
        const bool is_eq_sf_range = (_surface_range == rhs._surface_range);
        const bool is_eq_pt_range = (_portal_range == rhs._portal_range);
        const bool is_eq_sf_trf_range = (_surface_trf_range == rhs._surface_trf_range);
        const bool is_eq_pt_trf_range = (_portal_trf_range == rhs._portal_trf_range);
        const bool is_eq_sf_finder = (_surfaces_finder_entry == rhs._surfaces_finder_entry);

        return (is_eq_bounds && is_eq_index && is_eq_sf_range && is_eq_pt_range && is_eq_sf_trf_range && is_eq_pt_trf_range && is_eq_sf_finder);
    }

    /** Equality operator of volumes, convenience function - const
     *
     * @param rhs is the volume to be compared with
     *
     * checks identity for all index ranges and bounds
     **/
    bool operator==(const volume &rhs) {
        const bool is_eq_bounds = (_bounds == rhs._bounds);
        const bool is_eq_index = (_index == rhs._index);
        const bool is_eq_sf_range = (_surface_range == rhs._surface_range);
        const bool is_eq_pt_range = (_portal_range == rhs._portal_range);
        const bool is_eq_sf_trf_range = (_surface_trf_range == rhs._surface_trf_range);
        const bool is_eq_pt_trf_range = (_portal_trf_range == rhs._portal_trf_range);
        const bool is_eq_sf_finder = (_surfaces_finder_entry == rhs._surfaces_finder_entry);

        return (is_eq_bounds && is_eq_index && is_eq_sf_range && is_eq_pt_range && is_eq_sf_trf_range && is_eq_pt_trf_range && is_eq_sf_finder);
     }

    /** @return the bounds - const access */
    inline const array_type<scalar, 6> &bounds() const { return _bounds; }

    /** @return the name */
    inline const std::string &name() const { return _name; }

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

    /** Set the index into the detector surface container
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

    /** @return range of portals - const access */
    // const auto &portal_range() const { return _portal_range; }

    /** @return range of surface transforms - const access */
    template <bool surface_range = true>
    inline const auto &trf_range() const {
        if constexpr (surface_range) {
            return _surface_trf_range;
        } else {
            return _portal_trf_range;
        }
    }

    /** Set the index into the detector transform store for portals
     *
     * @param range Portal transform index range
     */
    template <bool surface_range = true>
    inline void set_trf_range(dindex_range range) {
        if constexpr (surface_range) {
            update_range(_surface_trf_range, std::move(range));
        } else {
            update_range(_portal_trf_range, std::move(range));
        }
    }

    /** Print volume.
     *
     * @returns the volume description as a string
     */
    inline const std::string to_string() const {
        std::stringstream ss;
        ss << " - name: '" << name() << "'" << std::endl;

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
    /** Volume section: name */
    std::string _name = "unknown";

    /** Bounds section, default for r, z, phi */
    array_type<scalar, 6> _bounds = {0.,
                                     std::numeric_limits<scalar>::max(),
                                     -std::numeric_limits<scalar>::max(),
                                     std::numeric_limits<scalar>::max(),
                                     -M_PI,
                                     M_PI};

    /** Volume index */
    dindex _index = dindex_invalid;

    /** Transform ranges in the detector transform store.*/
    dindex_range _surface_trf_range = {dindex_invalid, dindex_invalid};
    dindex_range _portal_trf_range = {dindex_invalid, dindex_invalid};

    /** Index ranges in the detector surface/portal containers.*/
    dindex_range _surface_range = {dindex_invalid, dindex_invalid};
    dindex_range _portal_range = {dindex_invalid, dindex_invalid};

    /** Index into the surface finder container */
    dindex _surfaces_finder_entry = dindex_invalid;

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
};

}  // namespace detray