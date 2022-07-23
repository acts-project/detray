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

/// @brief Base class for all masks.
///
/// Masks are lightweight types that define any special behaviour of different
/// surfaces geometries: They give the surface class its 'shape'.
/// A mask defines a local coordinate system on a surface geometry (e.g. a
/// plane surface or a cylinder surface) and defines the surface's extent in
/// that geometry. Masks are consequently used to make the inside/outside
/// decision during an intersection operation.
/// Masks can be (shared) boundary surfaces for detector volumes and therefore
/// hold the links to adjacent volumes that are needed during navigation. For
/// sensitive surfaces those links point to the mother volume.
///
/// @note Might contain multiple masks in the future.
///
/// @tparam intersector_t the intersection implementation for the surface
///                       geometry the mask belongs to. Can also be different,
///                       if the intersection is not done with a straight line.
/// @tparam locat_t the local coordinate system in which the surface extent is
///                 defined.
/// @tparam links_t the type of links that are used to tie detector volumes
///                 together. Might be multiindex in the future.
/// @tparam kDIM the number of values needed to define a surfaces extent.
template <typename intersector_t, typename local_t, typename links_t,
          template <typename, std::size_t> class array_t, std::size_t kDIM>
class mask_base {
    public:
    // Make template parameters accessible
    using intersector_type = intersector_t;
    using local_type = local_t;
    using links_type = links_t;

    template <typename T, std::size_t I>
    using array_type = array_t<T, I>;
    using mask_values = array_type<scalar, kDIM>;

    /// Default constructor
    mask_base() = default;

    DETRAY_HOST_DEVICE
    mask_base(const mask_values& val, const links_type& links)
        : _values(val), _volume_link(links) {}

    /// Equality operator from an array, convenience function
    ///
    /// @param rhs is the mask to be compared with
    ///
    /// @returns true if equal
    DETRAY_HOST_DEVICE
    bool operator==(const mask_values& rhs) { return (_values == rhs); }

    /// Equality operator
    ///
    /// @param rhs is the mask to be compared with
    ///
    /// @returns true if equal
    DETRAY_HOST_DEVICE
    bool operator==(const mask_base& rhs) {
        return (_values == rhs._values && _volume_link == rhs._volume_link);
    }

    /// Access operator - non-const
    /// @returns the reference to the requested mask value
    DETRAY_HOST_DEVICE
    scalar& operator[](std::size_t value_index) { return _values[value_index]; }

    /// Access operator - const
    /// @returns the reference to the requested mask value
    DETRAY_HOST_DEVICE
    scalar operator[](std::size_t value_index) const {
        return _values[value_index];
    }

    /// @return the mask values
    DETRAY_HOST_DEVICE
    const mask_values& values() const { return _values; }

    /// @return an associated intersector functor
    DETRAY_HOST_DEVICE
    intersector_type intersector() const { return intersector_type{}; };

    /// @return the local frame type
    DETRAY_HOST_DEVICE
    constexpr local_type local() const { return local_type{}; }

    /// @return the volume link - const reference
    DETRAY_HOST_DEVICE
    auto volume_link() const { return _volume_link; }

    /// @return the volume link - non-const access
    DETRAY_HOST_DEVICE
    auto volume_link() { return _volume_link; }

    protected:
    mask_values _values;
    links_type _volume_link;
};

}  // namespace detray