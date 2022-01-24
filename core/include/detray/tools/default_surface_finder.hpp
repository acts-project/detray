
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/** A dummy neighborhood finder that assumes all objects in a volume are stored
 *  contiguously.
 *
 * @tparam volume_t the volume that calls this surface finder
 * @tparam range_t the range indexing for objects
 */
template <typename volume_t, typename range_t =,
          template <typename...> class vector_type = dvector>
struct dummy_surface_finder {
    bool _sort = true;

    using point2 = __plugin::point2<detray::scalar>;

    /** Call operator for the local object search: Returns all objects.
     *
     * @param p2 the local 2d point. Not used
     *
     * @note the surface range
     **/
    auto operator()(const point2 & /*p2*/, const volume_t &vol) const {
        return vol.range();
    }
};

/** A brute force neighborhood finder that returns the indices for all objects
 *  in a volume.
 *
 * @tparam volume_t the volume that calls this surface finder
 */
template <typename volume_t, template <typename...> class vector_type = dvector>
struct brute_force_surface_finder {
    bool _sort = true;

    using point2 = __plugin::point2<detray::scalar>;

    /** Constructor from from volume
     *
     * @param vol the volume this surface finder works in
     **/
    brute_force_surface_finder(volume_t &vol) : {
        const auto &ranges = vol.ranges();
        for (const auto &range : ranges) {
            for (const auto &sf_idx : range) {
                _surface_indices.push_back(sf_idx);
            }
        }
    }

    /** Call operator for the indices of all objects inside the volume
     *
     * @param p2 the local 2d point for the search
     *
     * @return returns the vector of all indices
     **/
    auto operator()(const point2 &p2, const volume_t & /*vol*/) const {
        return _surface_indices;
    }

    vector_type<dindex> _surface_indices = {};
};

}  // namespace detray