/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::intersection {

namespace detail {

/// Filter the possible intersector types by the underlying geometry of the
/// detector surfaces (deduced from the mask shapes)
template <template <typename, typename, bool> class intersector_t,
          bool do_debug = !intersection::contains_pos>
struct select_intersector {
    // For a planar object, always assume a cartesian frame (e.g. rectangle
    // shape)
    template <typename mask_t>
    using type = std::conditional_t<
        concepts::planar_object<mask_t>,
        intersector_t<rectangle2D, typename mask_t::algebra_type, do_debug>,
        intersector_t<typename mask_t::shape, typename mask_t::algebra_type,
                      do_debug>>;
};

template <std::size_t N>
struct initialize {

    template <typename mask_group_t, typename mask_range_t,
              typename is_container_t, typename result_t, typename traj_t,
              typename surface_t, detray::concepts::transform3D transform_t,
              detray::concepts::scalar scalar_t>
    DETRAY_HOST_DEVICE DETRAY_INLINE constexpr void operator()(
        const mask_group_t &mask_group, const mask_range_t &mask_range,
        is_container_t &is_container, const result_t &result,
        const traj_t &traj, const surface_t &sf_desc, const transform_t &ctf,
        const intersection::config &cfg,
        const scalar_t external_mask_tolerance = 0.f) const {

        using intersection_t = typename is_container_t::value_type;

        // Resolve the masks that belong to the surface
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {

            intersection_t is{};

            // Build the resulting intersecion(s) from the intersection point
            if constexpr (N > 1) {
                std::uint8_t n_found{0u};

                for (std::size_t i = 0u; i < N; ++i) {
                    resolve_mask(is, traj, result[i], sf_desc, mask, ctf, cfg,
                                 external_mask_tolerance);

                    if (is.is_probably_inside()) {
                        insert_sorted(is, is_container);
                        ++n_found;
                    }
                    if (n_found == N) {
                        return;
                    }
                }
            } else {
                resolve_mask(is, traj, result, sf_desc, mask, ctf, cfg,
                             external_mask_tolerance);

                if (is.is_probably_inside()) {
                    insert_sorted(is, is_container);
                    return;
                }
            }
        }
    }

    template <typename intersection_t>
    DETRAY_HOST_DEVICE constexpr void insert_sorted(
        const intersection_t &sfi,
        std::vector<intersection_t> &intersections) const {

        auto itr_pos = detray::upper_bound(intersections.cbegin(),
                                           intersections.cend(), sfi);

        intersections.insert(itr_pos, sfi);
    }

    /// Specialization for the navigation state cache
    template <typename nav_state_t>
    DETRAY_HOST_DEVICE constexpr void insert_sorted(
        const typename nav_state_t::value_type &sfi,
        nav_state_t &intersections) const {

        auto itr_pos{intersections.cbegin()};

        // For just two candidates int the cache, the navigation state keeps
        // the first as the previouly visited candidate -> no sorting needed
        if constexpr (nav_state_t::capacity() > 2u) {
            itr_pos = detray::upper_bound(intersections.cbegin(),
                                          intersections.cend(), sfi);
        }

        intersections.insert(itr_pos, sfi);
    }
};

struct visitor {
    /// Operator function to calculate surface intersections
    ///
    /// @tparam intersector_t the intersector type corresponding to
    ///                       the surface mask shape id
    ///
    /// @param traj is the input trajectory
    /// @param surface is the input surface descriptor
    /// @param det is the detector data
    /// @param ctx is the geometry context
    /// @param cfg the config data for the intersector
    ///
    /// @return the number of valid intersections
    template <typename intersector_t, typename traj_t, typename is_container_t,
              typename detector_t>
    DETRAY_HOST_DEVICE DETRAY_INLINE constexpr auto operator()(
        const intersector_t & /*intersector*/, const traj_t &traj,
        const typename detector_t::surface_type &sf_desc,
        const typename detector_t::transform3_type &ctf,
        is_container_t &is_container, const detector_t &det,
        const intersection::config &cfg,
        const typename detector_t::scalar_type external_mask_tolerance =
            0.f) const {

        // Run the intersector
        constexpr const intersector_t intersector{};
        typename intersector_t::result_type result =
            intersector.point_of_intersection(traj, ctf,
                                              cfg.overstep_tolerance);

        // Check if any valid solutions were found
        if constexpr (intersector_t::n_solutions > 1) {
            bool found_any{false};
            for (const auto &ip : result) {
                if (ip.is_valid()) {
                    found_any = true;
                }
            }
            if (!found_any) [[unlikely]] {
                return;
            }
        } else {
            if (!result.is_valid()) [[unlikely]] {
                return;
            }
        }

        // Branch by mask
        det.mask_store()
            .template visit<detail::initialize<intersector_t::n_solutions>>(
                sf_desc.mask(), is_container, result, traj, sf_desc, ctf, cfg,
                external_mask_tolerance);
    }

    /// Specialization for the cylinder2D shape
    template <typename intersector_t, typename traj_t, typename is_container_t,
              typename detector_t>
        requires(
            std::same_as<typename intersector_t::frame_type,
                         cylindrical2D<typename detector_t::algebra_type>> &&
            detray::concepts::has_mask<detector_t, cylinder2D>)
    DETRAY_HOST_DEVICE DETRAY_INLINE constexpr auto operator()(
        const intersector_t & /*intersector*/, const traj_t &traj,
        const typename detector_t::surface_type &sf_desc,
        const typename detector_t::transform3_type &ctf,
        is_container_t &is_container, const detector_t &det,
        const intersection::config &cfg,
        const typename detector_t::scalar_type external_mask_tolerance =
            0.f) const {

        const auto &mask_coll =
            det.mask_store().template get<detector_t::masks::id::e_cylinder2>();

        const auto result = intersect_with_mask<intersector_t>(
            traj, ctf, sf_desc.mask().index(), mask_coll, cfg);

        // Check if any valid solutions were found
        bool found_any{false};
        for (const auto &ip : result) {
            if (ip.is_valid()) {
                found_any = true;
            }
        }
        if (!found_any) [[unlikely]] {
            return;
        }

        initialize<intersector_t::n_solutions>{}(
            mask_coll, sf_desc.mask().index(), is_container, result, traj,
            sf_desc, ctf, cfg, external_mask_tolerance);
    }

    /// Specialization for the concentric cylinder2D shape (cylinder portal)
    template <typename intersector_t, typename traj_t, typename is_container_t,
              typename detector_t>
        requires(
            std::same_as<
                typename intersector_t::frame_type,
                concentric_cylindrical2D<typename detector_t::algebra_type>> &&
            detray::concepts::has_mask<detector_t, concentric_cylinder2D>)
    DETRAY_HOST_DEVICE DETRAY_INLINE constexpr auto operator()(
        const intersector_t & /*intersector*/, const traj_t &traj,
        const typename detector_t::surface_type &sf_desc,
        const typename detector_t::transform3_type &ctf,
        is_container_t &is_container, const detector_t &det,
        const intersection::config &cfg,
        const typename detector_t::scalar_type external_mask_tolerance =
            0.f) const {

        const auto &mask_coll =
            det.mask_store()
                .template get<detector_t::masks::id::e_portal_cylinder2>();

        const auto result = intersect_with_mask<intersector_t>(
            traj, ctf, sf_desc.mask().index(), mask_coll, cfg);

        // Check if any valid solutions were found
        if (!result.is_valid()) [[unlikely]] {
            return;
        }

        initialize<intersector_t::n_solutions>{}(
            mask_coll, sf_desc.mask().index(), is_container, result, traj,
            sf_desc, ctf, cfg, external_mask_tolerance);
    }

    private:
    /// Operator function to calculate surface intersections with portal
    /// cylinders
    template <typename intersector_t, typename traj_t,
              detray::concepts::transform3D transform_t, typename mask_range_t,
              typename mask_coll_t>
    DETRAY_HOST_DEVICE DETRAY_INLINE constexpr auto intersect_with_mask(
        const traj_t &traj, const transform_t &ctf,
        const mask_range_t mask_range, const mask_coll_t &mask_coll,
        const intersection::config &cfg) const {

        // Get the mask index in this detector
        std::size_t mask_idx{detray::detail::invalid_value<std::size_t>()};
        if constexpr (detray::concepts::interval<mask_range_t>) {
            mask_idx = mask_range.lower();
        } else {
            mask_idx = mask_range;
        }

        assert(mask_idx < maks_coll.size());

        constexpr const intersector_t intersector{};
        return intersector.point_of_intersection(traj, ctf, mask_coll[mask_idx],
                                                 cfg.overstep_tolerance);
    }
};

template <template <typename, typename, bool> class intersector_t,
          bool do_debug = !intersection::contains_pos, typename traj_t,
          typename detector_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr auto calculate_intersection(
    const traj_t &traj, const typename detector_t::surface_type &sf_desc,
    const typename detector_t::transform3_type &ctf, const detector_t &det,
    const typename detector_t::geometry_context &ctx,
    const intersection::config &cfg) {
    // Filter the effective intersector types needed for the different
    // mask shapes and match them to the detector mask ids at compile-time
    using intersector_registry_t = types::mapped_registry<
        typename detector_t::masks,
        detail::select_intersector<intersector_t, do_debug>>;

    /*types::print<types::list<typename intersector_registry_t::type_list>>();

    for (auto i : intersector_registry_t::index_map()) {
        DETRAY_INFO_HOST(i);
    }

    types::print<types::list<typename detector_t::masks::type_list>>();*/

    static_assert(types::size<intersector_registry_t> <=
                  types::size<typename detector_t::masks>);

    // Visit only one of the intersectors according to the mask id
    return types::visit<intersector_registry_t, visitor>(
        sf_desc.mask().id(), traj, sf_desc, ctf, det, cfg);
}

}  // namespace detail

template <template <typename, typename, bool> class intersector_t,
          bool do_debug = !intersection::contains_pos, typename is_container_t,
          typename traj_t, typename detector_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr void intersect(
    const traj_t &traj, const typename detector_t::surface_type &sf_desc,
    is_container_t &is_container, const detector_t &det,
    const typename detector_t::geometry_context &ctx,
    const intersection::config &cfg,
    const typename detector_t::scalar_type external_mask_tolerance = 0.f) {

    const auto &ctf = det.transform_store().at(sf_desc.transform(), ctx);

    // Filter the effective intersector types needed for the different
    // mask shapes and match them to the detector mask ids at compile-time
    using intersector_registry_t = types::mapped_registry<
        typename detector_t::masks,
        detail::select_intersector<intersector_t, do_debug>>;

    static_assert(types::size<intersector_registry_t> <=
                  types::size<typename detector_t::masks>);

    // Visit only one of the intersectors according to the mask id
    types::visit<intersector_registry_t, detail::visitor>(
        sf_desc.mask().id(), traj, sf_desc, ctf, is_container, det, cfg,
        external_mask_tolerance);
}

/*template <template <typename, typename, bool> class intersector_t,
          bool do_debug = !intersection::contains_pos, typename intersection_t,
          typename traj_t, typename detector_t>
update(const traj_t &traj, intersection_t &is, const detector_t &det,
       const typename detector_t::geometry_context &ctx,
       const intersection::config &cfg,
       const typename detector_t::scalar_type external_mask_tolerance = 0.f) {

    // Branch by intersector
    const bool success =
        detail::calculate_intersection<intersector_t, do_debug>(
            traj, is.sf_desc, det, ctx, cfg);

    // Branch by mask
    if (success) {
        det.masks().template visit<update>();
    }
}*/

}  // namespace detray::intersection
