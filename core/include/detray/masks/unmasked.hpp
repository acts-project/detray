/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <string>

#include "detray/intersection/intersection.hpp"
#include "detray/intersection/planar_intersector.hpp"
#include "detray/intersection/unbound.hpp"

namespace detray {

template <typename intersector_t = planar_intersector,
          typename mask_local_t = unbound, typename mask_links_t = dindex>
struct unmasked {

    using mask_tolerance = bool;
    using local_type = mask_local_t;
    using links_type = mask_links_t;

    static constexpr mask_tolerance within_epsilon = true;

    links_type _links;

    /* Default constructor */
    unmasked() = default;

    /** Construction from voume link */
    DETRAY_HOST_DEVICE
    unmasked(links_type links) : _links(links) {}

    /** Mask operation
     *
     * @tparam point_type is the type of the point to be checked w.r.t. to
     * the mask bounds
     *
     * the parameters are ignored
     *
     * @return an intersection status e_inside / e_outside
     **/
    template <typename local_t>
    DETRAY_HOST_DEVICE inline intersection_status is_inside(
        const point2 & /*ignored*/,
        const mask_tolerance &t = within_epsilon) const {
        return t ? e_inside : e_outside;
    }

    /** Mask operation
     *
     * @tparam point_type is the type of the point to be checked w.r.t. to
     * the mask bounds
     *
     * the parameters are ignored
     *
     * @return true
     **/
    template <typename local_t>
    DETRAY_HOST_DEVICE inline bool is_inside(const point2 & /*ignored*/,
                                             scalar /*ignored*/) const {
        return true;
    }

    /** Mask operation
     *
     * @tparam point_type is the type of the point to be checked w.r.t. to
     * the mask bounds
     *
     * the parameters are ignored
     *
     * @return true
     **/
    template <typename local_t>
    DETRAY_HOST_DEVICE inline bool is_inside(const point2 & /*ignored*/,
                                             scalar /*ignored*/,
                                             scalar /*ignored*/) const {
        return true;
    }

    /** @return an associated intersector type */
    DETRAY_HOST_DEVICE
    intersector_t intersector() const { return intersector_t{}; };

    /** @return the local frame type */
    DETRAY_HOST_DEVICE
    constexpr local_type local() const { return local_type{}; }

    /** @return the links - const reference */
    DETRAY_HOST_DEVICE
    const links_type &links() const { return _links; }

    /** @return the links - non-const access */
    DETRAY_HOST_DEVICE
    links_type &links() { return _links; }

    /** @return the volume link - const reference */
    DETRAY_HOST_DEVICE
    dindex volume_link() const { return detail::get<0>(_links); }

    /** @return the volume link - non-const access */
    DETRAY_HOST_DEVICE
    dindex volume_link() { return detail::get<0>(_links); }

    /** @return the surface finder link - const reference */
    DETRAY_HOST_DEVICE
    dindex finder_link() const { return detail::get<1>(_links); }

    /** @return the surface finder link - non-const access */
    DETRAY_HOST_DEVICE
    dindex finder_link() { return detail::get<1>(_links); }

    /** Transform to a string for output debugging */
    std::string to_string() const { return "unmasked"; }
};

}  // namespace detray
