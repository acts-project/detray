/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/masks/masks.hpp"

// System include(s)
#include <cassert>
#include <limits>
#include <vector>

namespace detray {

/// An axis aligned bounding box of a given @tparam shape_t
template <typename shape_t = cuboid3D<>, typename scalar_t = scalar>
class axis_aligned_bounding_box {
    public:
    /// Define geometric properties
    /// @{
    using shape = shape_t;
    using axes = typename shape_t::template axes<>;
    template <typename algebra_t>
    using local_frame = typename shape_t::template local_frame_type<algebra_t>;

    static constexpr std::size_t Dim{axes::Dim};
    /// @}

    /// Default constructor builds an infinite box
    constexpr axis_aligned_bounding_box() = default;

    /// Constructor from single mask boundary values
    template <typename... Args>
    DETRAY_HOST_DEVICE explicit constexpr axis_aligned_bounding_box(
        unsigned int box_id, Args&&... args)
        : m_mask(box_id, std::forward<Args>(args)...) {}

    /// Constructor from a surface mask
    template <typename link_t, typename algebra_t,
              template <typename, std::size_t> class array_t>
    DETRAY_HOST_DEVICE explicit constexpr axis_aligned_bounding_box(
        const mask<shape, link_t, algebra_t, array_t>& mask,
        unsigned int box_id, const scalar_t envelope = 1.01f) {
        // Make sure the box is actually 'bounding'
        assert(envelope > (1.f + std::numeric_limits<scalar_t>::epsilon()));
        // Set this instances id
        m_mask.volume_link() = box_id;
        // Apply the envelope to the boundary values, which can be negative
        for (unsigned int i{0u}; i < shape::boundaries::e_size; ++i) {
            m_mask[i] = mask[i] * envelope;
        }
    }

    /// Constructor from mask boundary vector
    DETRAY_HOST axis_aligned_bounding_box(const std::vector<scalar_t>& values,
                                          unsigned int box_id)
        : m_mask(values, box_id) {
        assert(values.size() == shape::boundaries::e_size &&
               " Given number of boundaries does not match mask shape.");
    }

    /// Construct a bounding box around a set of boxes
    constexpr axis_aligned_bounding_box(
        const std::vector<axis_aligned_bounding_box>& aabbs) {}

    /// @returns the bounds of the box, depending on its shape
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> unsigned int { return m_mask.volume_link(); }

    /// @returns the bounds of the box, depending on its shape
    DETRAY_HOST_DEVICE
    constexpr auto bounds() const -> const mask<shape, unsigned int>& {
        return m_mask;
    }

    /// Checks whether a point lies inside the box. The point has to be defined
    /// in the coordinate frame that is spanned by the box axes.
    DETRAY_HOST_DEVICE
    template <typename point_t>
    constexpr auto is_inside(
        const point_t& loc_p,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const
        -> intersection::status {
        return m_mask.is_inside(loc_p, t);
    }

    /// Intersect the box with a ray
    /*DETRAY_HOST_DEVICE
    template<typename algebra_t>
    constexpr auto intersect(
        const detail::ray<algebra_t> &ray,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const
        -> intersection::status {
        static_assert(std::is_same_v<shape, cuboid3D<>>,
            "aabbs are only implemented in cuboid shape for now");
        return m_mask.intersector()(ray, m_mask, t);
    }*/

    /*DETRAY_HOST_DEVICE
    template<typename algebra_t>
    constexpr auto intersect(
        const axis_aligned_bounding_box &aabb,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const
        -> intersection::status {
        return m_mask.is_overlap(aabb.bounds());
    }*/

    /*DETRAY_HOST_DEVICE
    template<typename algebra_t>
    constexpr auto intersect(
        const detail::frustum<algebra_t> frustum,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const
        -> intersection::status {
        ....
    }*/

    private:
    /// Keeps the box boundary values and id
    mask<shape, unsigned int> m_mask;
};

}  // namespace detray