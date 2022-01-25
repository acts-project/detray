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
          template <typename, unsigned int> class array_t = darray>
class volume {

    public:
    // The type of objects a volume can contain
    using objects = object_registry_t;
    using obj_link_t = typename objects::link_type;
    using obj_range_t = typename objects::range_type;
    using sf_finders = sf_finder_registry_t;
    using sf_finder_link_t = typename sf_finders::link_type;
    // In case the geometry needs to be printed
    using name_map = std::map<dindex, std::string>;
    // used for sfinae
    using volume_def = volume<objects, sf_finders, array_t>;

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

    /** set surface finder link associated with volume */
    void set_surfaces_finder(sf_finder_link_t link) { _sf_finder = link; }

    /** get surface finder link associated with volume */
    const auto sf_finder_link() const { return _sf_finder; }

    /** get the type of surface finder associated with volume */
    const typename sf_finders::id sf_finder_type() const {
        return detail::get<0>(_sf_finder);
    }

    /** get the index of the surface finder associated with volume */
    const auto sf_finder_index() const { return detail::get<1>(_sf_finder); }

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
    DETRAY_HOST inline void update_range(const obj_range_t &other) {
        auto &rg = detail::get<range_id>(_ranges);
        // Range not set yet - initialize
        constexpr obj_range_t empty{};
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
        constexpr auto index = objects::template get_id<object_t>();
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
    template <typename obj_range_t>
    DETRAY_HOST_DEVICE inline dindex n_in_range(obj_range_t &&rg) const {
        return detail::get<1>(rg) - detail::get<0>(rg);
    }

    /** Stub for the surface finder call.
     * @return the number of elements in a given range
     */
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
                _ranges == rhs._ranges && _sf_finder == rhs._sf_finder);
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
    array_t<obj_range_t, objects::n_types> _ranges = {};

    /** Links to a specific surface finder for this volume */
    sf_finder_link_t _sf_finder = {};
};

}  // namespace detray