/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/cylindrical3D.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/masks/masks.hpp"

// System include(s)
#include <cassert>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace detray {

/// An axis aligned bounding box of a given @tparam shape_t
template <typename shape_t, typename T = detray::scalar,
          template <typename> class algebra_t = ALGEBRA_PLUGIN>
class axis_aligned_bounding_volume {
    using transform3D = dtransform3D<algebra_t<T>>;
    using point3D = dpoint3D<algebra_t<T>>;

    public:
    /// Define geometric properties
    /// @{
    using shape = shape_t;
    using boundaries = typename shape_t::boundaries;
    using axes = typename shape_t::template axes<>;
    using local_frame =
        typename shape_t::template local_frame_type<algebra_t<T>>;

    static constexpr std::size_t Dim{axes::Dim};
    /// @}

    /// Default constructor builds an infinite box
    constexpr axis_aligned_bounding_volume() = default;

    /// Constructor from mask boundary values
    template <typename... Args>
    DETRAY_HOST_DEVICE explicit constexpr axis_aligned_bounding_volume(
        std::size_t box_id, Args&&... args)
        : m_mask(box_id, std::forward<Args>(args)...) {}

    /// Construct around an arbitrary surface @param mask
    template <
        typename mask_t, typename S = shape_t,
        typename std::enable_if_t<std::is_same_v<S, cuboid3D<>>, bool> = true>
    DETRAY_HOST_DEVICE constexpr axis_aligned_bounding_volume(
        const mask_t& mask, std::size_t box_id, const T envelope)
        : m_mask{mask.local_min_bounds(envelope).values(), box_id} {
        // Make sure the box is actually 'bounding'
        assert(envelope >= std::numeric_limits<T>::epsilon());
    }

    /// Construct from mask boundary vector
    DETRAY_HOST axis_aligned_bounding_volume(const std::vector<T>& values,
                                             std::size_t box_id)
        : m_mask(values, box_id) {
        assert(values.size() == shape::boundaries::e_size &&
               " Given number of boundaries does not match mask shape.");
    }

    /// Construct a bounding box around a set of boxes
    /// @note the given bounding volumes need to be defnined in the same
    /// local coordinate system!
    template <typename other_shape_t, typename other_scalar_t,
              typename std::enable_if_t<
                  std::is_same_v<
                      typename shape::template local_frame_type<algebra_t<T>>,
                      typename other_shape_t::template local_frame_type<
                          algebra_t<T>>>,
                  bool> = true>
    DETRAY_HOST constexpr axis_aligned_bounding_volume(
        const std::vector<
            axis_aligned_bounding_volume<other_shape_t, other_scalar_t>>& aabbs,
        std::size_t box_id, const T env) {

        // Find min/max extent of the local aabb in local coordinates
        constexpr T inf{std::numeric_limits<T>::infinity()};
        T min_x{inf}, min_y{inf}, min_z{inf}, max_x{-inf}, max_y{-inf},
            max_z{-inf};
        for (const auto& vol : aabbs) {
            const auto min_point = vol.template loc_min();
            const auto max_point = vol.template loc_max();

            // Check every coordinate of the points
            min_x = min_point[0] < min_x ? min_point[0] : min_x;
            min_y = min_point[1] < min_y ? min_point[1] : min_y;

            max_x = max_point[0] > max_x ? max_point[0] : max_x;
            max_y = max_point[1] > max_y ? max_point[1] : max_y;

            if (min_point.size() > 1) {
                min_z = min_point[2] < min_z ? min_point[2] : min_z;
                max_z = max_point[2] > max_z ? max_point[2] : max_z;
            }
        }
        m_mask = mask<shape, std::size_t>{box_id,      min_x - env, min_y - env,
                                          min_z - env, max_x + env, max_y + env,
                                          max_z + env};
    }

    /// Subscript operator @returns a single box boundary.
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const std::size_t i) const -> T {
        return m_mask[i];
    }

    /// @returns the bounds of the box, depending on its shape
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> std::size_t { return m_mask.volume_link(); }

    /// @returns the bounds of the box, depending on its shape
    DETRAY_HOST_DEVICE
    constexpr auto bounds() const -> const mask<shape, std::size_t>& {
        return m_mask;
    }

    /// @returns the minimum bounds of the volume in local coordinates
    DETRAY_HOST_DEVICE constexpr auto loc_min() const -> point3D {
        if constexpr (std::is_same_v<shape, cuboid3D<>>) {
            return {m_mask[cuboid3D<>::e_min_x], m_mask[cuboid3D<>::e_min_y],
                    m_mask[cuboid3D<>::e_min_z]};
        } else if constexpr (std::is_same_v<shape, cylinder3D>) {
            return {-m_mask[cylinder3D::e_max_r], m_mask[cylinder3D::e_min_phi],
                    m_mask[cylinder3D::e_min_z]};
        }

        // If the volume shape is not supported, return universal minimum
        assert(false);
        constexpr T inf{std::numeric_limits<T>::infinity()};
        return {-inf, -inf, -inf};
    }

    /// @returns the maximum bounds of the volume in local coordinates
    DETRAY_HOST_DEVICE constexpr auto loc_max() const -> point3D {
        if constexpr (std::is_same_v<shape, cuboid3D<>>) {
            return {m_mask[cuboid3D<>::e_max_x], m_mask[cuboid3D<>::e_max_y],
                    m_mask[cuboid3D<>::e_max_z]};
        } else if constexpr (std::is_same_v<shape, cylinder3D>) {
            return {m_mask[cylinder3D::e_max_r], m_mask[cylinder3D::e_max_phi],
                    m_mask[cylinder3D::e_max_z]};
        }

        // If the volume shape is not supported, return universal minimum
        // (or compilation error for 2D point)
        assert(false);
        constexpr T inf{std::numeric_limits<T>::infinity()};
        return {inf, inf, inf};
    }
    /// @returns the minimum bounds of the volume in global cartesian
    /// coordinates
    DETRAY_HOST_DEVICE constexpr auto glob_min(const transform3D& trf) const
        -> point3D {

        if constexpr (std::is_same_v<shape, cuboid3D<>>) {
            return trf.point_to_global(loc_min());
        } else if constexpr (std::is_same_v<shape, cylinder3D>) {
            return cylindrical3D<algebra_t<T>>{}.local_to_global(trf, m_mask,
                                                                 loc_min());
        }

        // If the volume shape is not supported, return universal minimum
        // (or compilation error for 2D point)
        assert(false);
        constexpr T inf{std::numeric_limits<T>::infinity()};
        return {-inf, -inf, -inf};
    }

    /// @returns the maximum bounds of the volume in global cartesian
    /// coordinates
    DETRAY_HOST_DEVICE constexpr auto glob_max(const transform3D& trf) const
        -> point3D {

        if constexpr (std::is_same_v<shape, cuboid3D<>>) {
            return trf.point_to_global(loc_max());
        } else if constexpr (std::is_same_v<shape, cylinder3D>) {
            return cylindrical3D<algebra_t<T>>{}.local_to_global(trf, m_mask,
                                                                 loc_max());
        }

        // If the volume shape is not supported, return universal minimum
        assert(false);
        constexpr T inf{std::numeric_limits<T>::infinity()};
        return {inf, inf, inf};
    }

    /// @returns the geometric center position in global cartesian system
    DETRAY_HOST_DEVICE constexpr auto center() const -> point3D {

        const T center_x{
            0.5f * (m_mask[cuboid3D<>::e_max_x] + m_mask[cuboid3D<>::e_min_x])};
        const T center_y{
            0.5f * (m_mask[cuboid3D<>::e_max_y] + m_mask[cuboid3D<>::e_min_y])};
        const T center_z{
            0.5f * (m_mask[cuboid3D<>::e_max_z] + m_mask[cuboid3D<>::e_min_z])};

        return {std::isinf(center_x) ? 0.f : center_x,
                std::isinf(center_y) ? 0.f : center_y,
                std::isinf(center_z) ? 0.f : center_z};
    }

    /// @brief Lower and upper point for minimum axis aligned bounding box of
    /// cuboid shape.
    ///
    /// Computes the min and max vertices in global coordinates from the local
    /// minimum aabb. The global aabb is not necessarily an minimum aabb.
    ///
    /// @param trf affine transformation
    ///
    /// @returns a new, transformed aabb.
    template <
        typename S = shape_t,
        typename std::enable_if_t<std::is_same_v<S, cuboid3D<>>, bool> = true>
    DETRAY_HOST_DEVICE auto transform(const transform3D& trf) const
        -> axis_aligned_bounding_volume {

        const T scalor_x{
            (m_mask[cuboid3D<>::e_max_x] - m_mask[cuboid3D<>::e_min_x])};
        const T scalor_y{
            (m_mask[cuboid3D<>::e_max_y] - m_mask[cuboid3D<>::e_min_y])};
        const T scalor_z{
            (m_mask[cuboid3D<>::e_max_z] - m_mask[cuboid3D<>::e_min_z])};

        // Cannot handle 'inf' propagation through the calculation for now
        if (std::isinf(scalor_x) or std::isinf(scalor_y) or
            std::isinf(scalor_z)) {
            // If the box was infinite to begin with, it stays that way
            assert(std::isinf(scalor_x) and std::isinf(scalor_y) and
                   std::isinf(scalor_z));

            return *this;
        }

        // The new axis vectors, scaled to the aabb dimensions
        // (e.g. max_x - min_x)
        const point3D new_box_x = scalor_x * trf.x();
        const point3D new_box_y = scalor_y * trf.y();
        const point3D new_box_z = scalor_z * trf.z();

        // Transform the old min and max points to the global frame and
        // construct all corner points of the local aabb in global coordinates
        std::array<point3D, 8> glob_c_points;
        glob_c_points[0] = glob_min(trf);
        glob_c_points[1] = glob_max(trf);
        glob_c_points[2] = glob_c_points[0] + new_box_x;
        glob_c_points[3] = glob_c_points[0] + new_box_y;
        glob_c_points[4] = glob_c_points[0] + new_box_z;
        glob_c_points[5] = glob_c_points[1] - new_box_x;
        glob_c_points[6] = glob_c_points[1] - new_box_y;
        glob_c_points[7] = glob_c_points[1] - new_box_z;

        // Find min/max extent of the local aabb in global coordinates
        constexpr T inf{std::numeric_limits<T>::infinity()};
        T min_x{inf}, min_y{inf}, min_z{inf}, max_x{-inf}, max_y{-inf},
            max_z{-inf};
        for (const point3D& p : glob_c_points) {
            // Check every coordinate of the point
            min_x = p[0] < min_x ? p[0] : min_x;
            max_x = p[0] > max_x ? p[0] : max_x;
            min_y = p[1] < min_y ? p[1] : min_y;
            max_y = p[1] > max_y ? p[1] : max_y;
            min_z = p[2] < min_z ? p[2] : min_z;
            max_z = p[2] > max_z ? p[2] : max_z;
        }

        // Construct the transformed aabb
        return axis_aligned_bounding_volume{
            m_mask.volume_link(), min_x, min_y, min_z, max_x, max_y, max_z};
    }

    /// Checks whether a point lies inside the box. The point has to be defined
    /// in the coordinate frame that is spanned by the box axes.
    DETRAY_HOST_DEVICE
    template <typename point_t>
    constexpr auto is_inside(
        const point_t& loc_p,
        const T t = std::numeric_limits<T>::epsilon()) const {
        return m_mask.is_inside(loc_p, t);
    }

    /// Intersect the box with a ray
    DETRAY_HOST_DEVICE
    constexpr bool intersect(
        const detail::ray<transform3D>& ray,
        const T t = std::numeric_limits<T>::epsilon()) const {
        static_assert(std::is_same_v<shape, cuboid3D<>>,
                      "aabbs are only implemented in cuboid shape for now");
        return m_mask
            .template intersector<intersection2D<bool, T, algebra_t>>()(
                ray, m_mask, t);
    }

    /// @TODO: Overlapping aabbs
    /*DETRAY_HOST_DEVICE
    constexpr auto intersect(
        const axis_aligned_bounding_volume &aabb,
        const T t = std::numeric_limits<T>::epsilon()) const
        {
        return m_mask.is_overlap(aabb.bounds());
    }*/

    /// @TODO: Frustum intersection
    /*DETRAY_HOST_DEVICE
    constexpr auto intersect(
        const detail::frustum<T, algebra_t> frustum,
        const T t = std::numeric_limits<T>::epsilon()) const
        {
        ....
    }*/

    private:
    /// Keeps the box boundary values and id
    mask<shape, std::size_t, algebra_t<T>> m_mask;
};

template <typename shape_t, typename T = scalar,
          template <typename> class algebra_t = ALGEBRA_PLUGIN>
DETRAY_HOST std::ostream& operator<<(
    std::ostream& os,
    const axis_aligned_bounding_volume<shape_t, T, algebra_t>& aabb) {
    return os << aabb.bounds().to_string();
}

}  // namespace detray
