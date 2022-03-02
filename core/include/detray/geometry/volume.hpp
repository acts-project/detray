/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

/** Volume class that acts as a logical container in the geometry for geometry
 *  objects, such as module surfaces or portals. The information is kept as
 *  index ranges into larger containers that are owned by the detector
 *  implementation.
 *
 * @tparam array_t the type of the internal array, must have STL semantics
 */
template <typename object_registry_t, typename range_t = dindex_range,
          template <typename, std::size_t> class array_t = darray>
class volume {

    public:
    // The type of objects a volume can contain
    using objects = object_registry_t;
    // In case the geometry needs to be printed
    using name_map = std::map<dindex, std::string>;
    // used for sfinae
    using volume_def = volume<object_registry_t, range_t, array_t>;

    enum grid_type : unsigned int {
        e_no_grid = 0,
        e_z_phi_grid = 1,  // barrel
        e_r_phi_grid = 2,  // endcap
    };

    /** Default constructor**/
    volume() = default;

    /** Contructor with bounds
     * @param bounds of the volume
     */
    volume(const array_t<scalar, 6> &bounds) : _bounds(bounds) {}

    /** @return the bounds - const access */
    DETRAY_HOST_DEVICE
    inline const array_t<scalar, 6> &bounds() const { return _bounds; }

    /** @return the name (add an offset for the detector name)*/
    DETRAY_HOST_DEVICE
    inline const std::string &name(const name_map &names) const {
        return names.at(_index + 1);
    }

    /** @return the index */
    DETRAY_HOST_DEVICE
    inline dindex index() const { return _index; }

    /** @param index the index */
    DETRAY_HOST
    inline void set_index(const dindex index) { _index = index; }

    /** @return the entry into the local surface finders */
    DETRAY_HOST_DEVICE
    inline dindex surfaces_finder_entry() const {
        return _surfaces_finder_entry;
    }

    /** @param entry the entry into the local surface finders */
    DETRAY_HOST
    inline void set_surfaces_finder(const dindex entry) {
        _surfaces_finder_entry = entry;
    }

    /** @return if the volume is empty or not */
    DETRAY_HOST_DEVICE inline bool empty() const {
        return n_objects<objects::e_surface>() == 0;
    }

    /** @return the number of surfaces in the volume */
    template <typename objects::id range_id = objects::e_surface>
    DETRAY_HOST_DEVICE inline auto n_objects() const {
        return n_in_range(range<range_id>());
    }

    /** Set or update the index into a geometry container identified by the
     *  range_id.
     *
     * @param other Surface index range
     */
    template <typename objects::id range_id = objects::e_surface>
    DETRAY_HOST inline void update_range(const range_t &other) {
        auto &rg = detail::get<range_id>(_ranges);
        // Range not set yet - initialize
        constexpr range_t empty{};
        if (rg == empty) {
            rg = other;
        } else {
            // Update
            assert(detail::get<1>(rg) == detail::get<0>(other));
            detail::get<1>(rg) = detail::get<1>(other);
        }
    }

    /** @return range of surfaces by surface type - const access */
    template <typename object_t>
    DETRAY_HOST_DEVICE inline constexpr const auto &range() const {
        constexpr auto index = objects::template get_index<object_t>::value;
        return detail::get<index>(_ranges);
    }

    /** @return range of surfaces- const access */
    template <typename objects::id range_id = objects::e_surface>
    DETRAY_HOST_DEVICE inline constexpr const auto &range() const {
        return detail::get<range_id>(_ranges);
    }

    /** @return _ranges */
    DETRAY_HOST_DEVICE inline constexpr const auto &ranges() const {
        return _ranges;
    }

    /** @return the number of elements in a given range */
    template <typename range_type>
    DETRAY_HOST_DEVICE inline dindex n_in_range(range_type &&rg) const {
        return detail::get<1>(rg) - detail::get<0>(rg);
    }

    /** set grid type associated with volume */
    void set_grid_type(grid_type val) { _grid_type = val; }

    /** get grid type associated with volume */
    auto get_grid_type() const { return _grid_type; }

    /** Equality operator
     *
     * @param rhs is the right hand side to be compared to
     */
    DETRAY_HOST_DEVICE
    bool operator==(const volume &rhs) const {
        return (_bounds == rhs._bounds && _index == rhs._index &&
                _ranges == rhs._ranges &&
                _surfaces_finder_entry == rhs._surfaces_finder_entry);
    }

    private:
    /** Bounds section, default for r, z, phi */
    array_t<scalar, 6> _bounds = {0.,
                                  std::numeric_limits<scalar>::max(),
                                  -std::numeric_limits<scalar>::max(),
                                  std::numeric_limits<scalar>::max(),
                                  -M_PI,
                                  M_PI};

    /** Volume index */
    dindex _index = dindex_invalid;

    /** Ranges in geometry containers for different objects types that belong
     * to this volume */
    array_t<range_t, objects::n_types> _ranges = {};

    /** Index into the surface finder container */
    dindex _surfaces_finder_entry = dindex_invalid;

    /** grid type associated with volume **/
    grid_type _grid_type = e_no_grid;
};

}  // namespace detray