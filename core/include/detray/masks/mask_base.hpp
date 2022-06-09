/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename intersector_t, typename local_t, typename links_t,
          template <typename, std::size_t> class array_t>
class mask_base {
    public:
    using intersector_type = intersector_t;
    using local_type = local_t;
    using links_type = links_t;

    template <typename T, std::size_t I>
    using array_type = array_t<T, I>;

    /// Return an associated intersector type
    DETRAY_HOST_DEVICE
    intersector_type intersector() const { return intersector_type{}; };

    /// Return the local frame type
    DETRAY_HOST_DEVICE
    constexpr local_type local() const { return local_type{}; }

    /// @return the links - const reference
    DETRAY_HOST_DEVICE
    const links_type &links() const { return _links; }

    /// @return the volume link - const reference
    DETRAY_HOST_DEVICE
    auto volume_link() const { return detail::get<0>(_links); }

    /// @return the volume link - non-const access
    DETRAY_HOST_DEVICE
    auto volume_link() { return detail::get<0>(_links); }

    /// @return the surface finder link - const reference
    DETRAY_HOST_DEVICE
    auto finder_link() const { return detail::get<1>(_links); }

    /// @return the surface finder link - non-const access
    DETRAY_HOST_DEVICE
    auto finder_link() { return detail::get<1>(_links); }

    protected:
    links_type _links;
};

}  // namespace detray