/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

/** Volume class that acts as a logical container in the geometry for geometry
 *  objects, such as module surfaces or portals. The information is kept as
 *  index ranges into larger containers that are owned by the detector
 *  implementation.
 *
 * @tparam array_type the type of the internal array, must have STL semantics
 */
template <typename object_registry_t, typename sf_finder_registry_t,
          typename range_t = dindex_range,
          template <typename, unsigned int> class array_t = darray,
          typename surface_finder_data_t = dindex>
class volume {

    public:
    // The type of objects a volume can contain
    using objects = object_registry_t;
    using sf_finders = sf_finder_registry_t;
    // In case the geometry needs to be printed
    using name_map = std::map<dindex, std::string>;
    // used for sfinae
    using volume_def =
        volume<objects, sf_finders, range_t, array_t, surface_finder_data_t>;

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
        return _surfaces_finder_data;
    }

    /** @param entry the entry into the local surface finders */
    DETRAY_HOST
    inline void set_surfaces_finder(const dindex entry) {
        _surfaces_finder_data = entry;
    }

    /** set grid type associated with volume */
    void set_grid_type(typename sf_finders::id val) { _grid_type = val; }

    /** get grid type associated with volume */
    const auto get_grid_type() const { return _grid_type; }

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
        auto &rg = std::get<range_id>(_ranges);
        // Range not set yet - initialize
        constexpr range_t empty{};
        if (rg == empty) {
            rg = other;
        } else {
            // Update
            assert(std::get<1>(rg) == std::get<0>(other));
            std::get<1>(rg) = std::get<1>(other);
        }
    }

    /** @return range of surfaces by surface type - const access */
    template <typename object_t>
    DETRAY_HOST_DEVICE inline constexpr const auto &range() const {
        constexpr auto index = objects::template get_id<object_t>();
        return std::get<index>(_ranges);
    }

    /** @return range of surfaces- const access */
    template <typename objects::id range_id = objects::e_surface>
    DETRAY_HOST_DEVICE inline constexpr const auto &range() const {
        return std::get<range_id>(_ranges);
    }

    /** @return _ranges */
    DETRAY_HOST_DEVICE inline constexpr const auto &ranges() const {
        return _ranges;
    }

    /** @return the number of elements in a given range */
    template <typename range_id>
    DETRAY_HOST_DEVICE inline dindex n_in_range(range_id &&rg) const {
        return std::get<1>(rg) - std::get<0>(rg);
    }

    /** @return the number of elements in a given range */
    template <typename point_t>
    DETRAY_HOST_DEVICE inline auto get_neighbors(point_t & /*pt*/) const {
        return range<objects::e_surface>();
    }

    /** Equality operator
     *
     * @param rhs is the right hand side to be compared to
     */
    DETRAY_HOST_DEVICE
    bool operator==(const volume &rhs) const {
        return (_bounds == rhs._bounds && _index == rhs._index &&
                _ranges == rhs._ranges &&
                _surfaces_finder_data == rhs._surfaces_finder_data);
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

    /** grid type associated with volume **/
    typename sf_finders::id _grid_type = sf_finders::id::e_unknown;

    /** Index into the surface finder container */
    surface_finder_data_t _surfaces_finder_data = dindex_invalid;
};

}  // namespace detray