/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/materials/detail/concepts.hpp"
#include "detray/materials/detail/material_accessor.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray::detail {

/// Functors to be used in the @c surface class
template <typename algebra_t>
struct surface_kernels {

    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    /// A functor to retrieve the masks shape name
    struct get_shape_name {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST inline std::string operator()(const mask_group_t&,
                                                  const index_t&) const {

            return std::string(mask_group_t::value_type::shape::name);
        }
    };

    /// A functor that checks if a local point @param loc_p is within the
    /// surface mask with tolerance @param tol
    struct is_inside {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline bool operator()(
            const mask_group_t& mask_group, const index_t& index,
            const point3_type& loc_p, const scalar_type tol) const {

            return mask_group[index].is_inside(loc_p, tol);
        }
    };

    /// A functor to run the mask self check. Puts error messages into @param os
    struct mask_self_check {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index,
            std::ostream& os) const {

            return mask_group.at(index).self_check(os);
        }
    };

    /// A functor to retrieve the masks volume link
    struct get_volume_link {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index) const {

            return mask_group[index].volume_link();
        }
    };

    /// A functor to retrieve a mask boundary, determined by @param i
    struct get_mask_value {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index,
            std::size_t i) const {

            return mask_group[index][i];
        }
    };

    /// A functor to retrieve the mask boundaries (host only)
    struct get_mask_values {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                           const index_t& index) const {

            std::vector<scalar_type> values{};
            for (const scalar_type v : mask_group[index].values()) {
                values.push_back(v);
            }

            return values;
        }
    };

    /// A functor to retrieve the material parameters
    struct get_material_params {
        template <typename mat_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mat_group_t& mat_group, const index_t& idx,
            const point2_type& loc_p) const {

            using material_t = typename mat_group_t::value_type;

            if constexpr (concepts::surface_material<material_t>) {
                return &(detail::material_accessor::get(mat_group, idx, loc_p)
                             .get_material());
            } else {
                using scalar_t = typename material_t::scalar_type;
                // Volume material (cannot be reached from a surface)
                return static_cast<const material<scalar_t>*>(nullptr);
            }
        }
    };

    /// A functor get the surface normal at a given local/bound position
    struct normal {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3_type operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3_type& trf3, const point2_type& bound) const {
            using mask_t = typename mask_group_t::value_type;

            return mask_t::get_local_frame().normal(trf3, bound,
                                                    mask_group[index]);
        }

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3_type operator()(
            const mask_group_t& /*mask_group*/, const index_t& /*index*/,
            const transform3_type& trf3, const point3_type& local) const {
            using mask_t = typename mask_group_t::value_type;

            return mask_t::get_local_frame().normal(trf3, local);
        }
    };

    /// A functor get the mask centroid in local cartesian coordinates
    struct centroid {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3_type operator()(
            const mask_group_t& mask_group, const index_t& index) const {

            return mask_group[index].centroid();
        }
    };

    /// A functor to perform global to local transformation
    struct global_to_local {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3_type operator()(
            const mask_group_t& /*mask_group*/, const index_t& /*index*/,
            const transform3_type& trf3, const point3_type& global,
            const vector3_type& dir) const {
            using mask_t = typename mask_group_t::value_type;

            return mask_t::to_local_frame(trf3, global, dir);
        }
    };

    /// A functor to perform local to global transformation
    struct local_to_global {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3_type operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3_type& trf3, const point2_type& bound,
            const vector3_type& dir) const {
            using mask_t = typename mask_group_t::value_type;

            return mask_t::get_local_frame().local_to_global(
                trf3, mask_group[index], bound, dir);
        }

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3_type operator()(
            const mask_group_t& /*mask_group*/, const index_t& /*index*/,
            const transform3_type& trf3, const point3_type& local,
            const vector3_type&) const {
            using mask_t = typename mask_group_t::value_type;

            return mask_t::to_global_frame(trf3, local);
        }
    };
};

}  // namespace detray::detail
