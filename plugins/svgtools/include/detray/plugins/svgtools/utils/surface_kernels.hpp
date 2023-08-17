/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/masks/masks.hpp"
#include "detray/plugins/svgtools/utils/volume_utils.hpp"

// System include(s)
#include <cassert>
#include <optional>

namespace detray::svgtools::utils {

/// @brief Functor to calculate the outermost radius a shape.
/// If the shape is not defined by a radius, then null option is returned.
struct outer_radius_getter {

    public:
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline std::optional<detray::scalar> operator()(
        const mask_group_t& mask_group, const index_t& index) const {
        const auto& m = mask_group[index];
        return outer_radius(m);
    }

    private:
    // The remaining shapes do not have an outer radius.
    template <typename mask_t>
    std::optional<detray::scalar> inline outer_radius(
        const mask_t& /*mask*/) const {
        return std::nullopt;
    }

    // Calculates the outer radius for rings.
    auto inline outer_radius(const detray::mask<detray::ring2D<>>& mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_outer_r]);
    }

    // Calculates the outer radius for annuluses.
    auto inline outer_radius(
        const detray::mask<detray::annulus2D<>>& mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_max_r]);
    }

    // Calculates the outer radius for cylinders (2D).
    template <bool kRadialCheck, template <typename> class intersector_t>
    auto inline outer_radius(
        const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>&
            mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_r]);
    }

    // Calculates the outer radius for cylinders (3D).
    auto inline outer_radius(
        const detray::mask<detray::cylinder3D>& mask) const {
        return std::optional<detray::scalar>(mask[mask.get_shape().e_max_r]);
    }
};

/// @brief Functor to obtain the volume link.
struct link_getter {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                       const index_t& index) const {
        const auto& m = mask_group[index];
        return m.volume_link();
    }
};

/// @brief Functor to calculate a suitable starting point for displaying the
/// link arrow.
struct link_start_getter {

    public:
    template <typename mask_group_t, typename index_t, typename transform_t>
    DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                       const index_t& index,
                                       const transform_t& transform) const {
        const auto& m = mask_group[index];
        return link_start(m, transform);
    }

    private:
    // Calculates the link starting location of the remaining shapes.
    template <typename mask_t, typename transform_t>
    auto inline link_start(const mask_t& mask,
                           const transform_t& transform) const {
        return mask.global_min_bounds_center(transform);
    }

    // Calculates the (optimal) link starting point for rings.
    template <typename transform_t>
    auto inline link_start(const detray::mask<detray::ring2D<>>& mask,
                           const transform_t& transform) const {
        using mask_t = typename detray::mask<detray::ring2D<>>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{(mask[shape_t::e_inner_r] + mask[shape_t::e_outer_r]) /
                         2.f};
        const scalar_t phi{detray::constant<scalar_t>::pi_2};
        const scalar_t z{0};

        // Polar coordinate system.
        const typename mask_t::local_frame_type frame{};

        return frame.local_to_global(transform, point3_t{r, phi, z});
    }

    // Calculates the (optimal) link starting point for annuluses.
    template <typename transform_t>
    auto inline link_start(const detray::mask<detray::annulus2D<>>& mask,
                           const transform_t& transform) const {
        using mask_t = typename detray::mask<detray::annulus2D<>>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{(mask[shape_t::e_min_r] + mask[shape_t::e_max_r]) / 2.f};
        const scalar_t phi{mask[shape_t::e_average_phi]};
        const scalar_t z{0};

        // Polar coordinate system.
        const typename mask_t::local_frame_type frame{};

        const auto true_center = mask.global_min_bounds_center(transform);
        const auto rel_point =
            frame.local_to_global(transform, point3_t{r, phi, z}) -
            transform.translation();
        return rel_point + true_center;
    }

    // Calculates the (optimal) link starting point for cylinders (2D).
    template <typename transform_t, bool kRadialCheck,
              template <typename> class intersector_t>
    auto inline link_start(
        const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>&
            mask,
        const transform_t& transform) const {
        using mask_t = typename detray::mask<
            detray::cylinder2D<kRadialCheck, intersector_t>>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{mask[shape_t::e_r]};
        const scalar_t phi{detray::constant<scalar_t>::pi_2};
        const scalar_t z{0};

        // Cylindrical coordinate system.
        const typename mask_t::local_frame_type frame{};

        const auto true_center = mask.global_min_bounds_center(transform);
        const auto rel_point =
            frame.local_to_global(transform, point3_t{r * phi, z, r}) -
            transform.translation();
        return rel_point + true_center;
    }

    // Calculates the (optimal) link starting point for cylinders (3D).
    template <typename transform_t>
    auto inline link_start(const detray::mask<detray::cylinder3D>& mask,
                           const transform_t& transform) const {
        using mask_t = typename detray::mask<detray::cylinder3D>;
        using shape_t = typename mask_t::shape;
        using point3_t = typename mask_t::point3_t;
        using scalar_t = typename mask_t::scalar_type;

        const scalar_t r{(mask[shape_t::e_min_r] + mask[shape_t::e_max_r]) / 2.f};
        const scalar_t phi{
            (mask[shape_t::e_max_phi] + mask[shape_t::e_max_phi]) / 2.f};
        const scalar_t z{0};

        // Cylindrical coordinate system.
        const typename mask_t::local_frame_type frame{};

        const auto true_center = mask.global_min_bounds_center(transform);
        const auto rel_point =
            frame.local_to_global(transform, point3_t{r * phi, z, r}) -
            transform.translation();
        return rel_point + true_center;
    }
};

/// @brief Functor to calculate a suitable end point for displaying the link
/// arrow.
struct link_end_getter {

    public:
    template <typename mask_group_t, typename index_t, typename detector_t,
              typename point3_t, typename vector3_t, typename scalar_t>
    DETRAY_HOST inline auto operator()(
        const mask_group_t& mask_group, const index_t& index,
        const detector_t& detector,
        const detray::detector_volume<detector_t>& volume,
        const point3_t& surface_point, const vector3_t& surface_normal,
        const scalar_t& link_length) const {
        const auto& m = mask_group[index];
        return link_dir(m, detector, volume, surface_point, surface_normal) *
                   link_length +
               surface_point;
    }

    private:
    /// @brief Calculates the direction of the link for remaining shapes.
    template <typename detector_t, typename mask_t, typename point3_t,
              typename vector3_t>
    inline auto link_dir(const mask_t& /*mask*/, const detector_t& /*detector*/,
                         const detray::detector_volume<detector_t>& volume,
                         const point3_t& surface_point,
                         const vector3_t& surface_normal) const {
        const auto dir = volume.center() - surface_point;
        const auto dot_prod = vector::dot(dir, surface_normal);
        // Should geometrically not happen with a local point 'surface_point'
        assert(dot_prod != 0.f);
        typename detector_t::scalar_type sgn =
            dot_prod > 0.f ? 1.f : (dot_prod < 0.f ? -1.f : 0);
        return sgn * surface_normal;
    }

    /// @brief Calculates the direction of the link for cylinders (2D)
    template <typename detector_t, bool kRadialCheck,
              template <typename> class intersector_t, typename point3_t,
              typename vector3_t>
    inline auto link_dir(
        const detray::mask<detray::cylinder2D<kRadialCheck, intersector_t>>&
            mask,
        const detector_t& detector,
        const detray::detector_volume<detector_t>& volume,
        const point3_t& /*surface_point*/,
        const vector3_t& surface_normal) const {
        for (const auto& desc :
             svgtools::utils::surface_lookup(detector, volume)) {
            const detray::surface surface{detector, desc};
            if (surface.is_portal()) {
                if (auto radius =
                        surface.template visit_mask<outer_radius_getter>()) {
                    if (*radius > mask[mask.get_shape().e_r]) {
                        return surface_normal;
                    }
                }
            }
        }
        return -1.f * surface_normal;
    }
};

}  // namespace detray::svgtools::utils