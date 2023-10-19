/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracks/tracks.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray::detail {

/// Functors to be used in the @c surface class
template <typename algebra_t>
struct surface_kernels {

    using transform3 = algebra_t;
    using point2 = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;
    using free_vector_type =
        typename free_track_parameters<algebra_t>::vector_type;
    using bound_vector_type =
        typename bound_track_parameters<algebra_t>::vector_type;
    using matrix_operator = typename algebra_t::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;

    /// A functor to retrieve the masks volume link
    struct get_volume_link {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index) const {

            return mask_group[index].volume_link();
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

    /// A functor get the surface normal at a given local/bound position
    struct normal {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point2& bound) const {

            const auto& m = mask_group[index];

            return m.local_frame().normal(trf3, bound, m);
        }

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point3& local) const {

            return mask_group[index].local_frame().normal(trf3, local);
        }
    };

    /// A functor to perform global to local bound transformation
    struct global_to_bound {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point2 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point3& global,
            const vector3& dir) const {

            const point3 local =
                mask_group[index].to_local_frame(trf3, global, dir);

            return {local[0], local[1]};
        }
    };

    /// A functor to perform global to local transformation
    struct global_to_local {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point3& global,
            const vector3& dir) const {

            return mask_group[index].to_local_frame(trf3, global, dir);
        }
    };

    /// A functor to perform local to global transformation
    struct local_to_global {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point2& bound,
            const vector3& dir) const {

            const auto& m = mask_group[index];

            return mask_group[index].local_frame().bound_local_to_global(
                trf3, m, bound, dir);
        }

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point3& local, const vector3&) const {

            return mask_group[index].to_global_frame(trf3, local);
        }
    };

    /// A functor to get from a free to a bound vector
    struct free_to_bound_vector {

        // Visitor to the detector mask store that is called on the mask
        // collection that contains the mask (shape) type of the surface
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline bound_vector_type operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const free_vector_type& free_vec) const {

            return mask_group[index].local_frame().free_to_bound_vector(
                trf3, free_vec);
        }
    };

    /// A functor to get from a bound to a free vector
    struct bound_to_free_vector {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline free_vector_type operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const bound_vector_type& bound_vec) const {

            const auto& m = mask_group[index];

            return m.local_frame().bound_to_free_vector(trf3, m, bound_vec);
        }
    };

    /// A functor to get the free-to-bound Jacobian
    struct free_to_bound_jacobian {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const free_vector_type& free_vec) const {

            return mask_group[index].local_frame().free_to_bound_jacobian(
                trf3, free_vec);
        }
    };

    /// A functor to get the bound-to-free Jacobian
    struct bound_to_free_jacobian {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const bound_vector_type& bound_vec) const {

            const auto& m = mask_group[index];

            return m.local_frame().bound_to_free_jacobian(trf3, m, bound_vec);
        }
    };

    /// A functor to get the path correction
    struct path_correction {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline free_matrix operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const vector3& pos, const vector3& dir,
            const vector3& dtds) const {

            return mask_group[index].local_frame().path_correction(pos, dir,
                                                                   dtds, trf3);
        }
    };

    /// A functor to get the local min bounds.
    struct local_min_bounds {

        template <typename mask_group_t, typename index_t, typename scalar_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index,
            const scalar_t env =
                std::numeric_limits<scalar_t>::epsilon()) const {

            return mask_group[index].local_min_bounds(env);
        }
    };

    /// A functor to get the vertices in local coordinates.
    struct local_vertices {

        template <typename mask_group_t, typename index_t, typename scalar_t>
        DETRAY_HOST_DEVICE inline auto operator()(
            const mask_group_t& mask_group, const index_t& index,
            const dindex n_seg) const {

            return mask_group[index].vertices(n_seg);
        }
    };
};

}  // namespace detray::detail
