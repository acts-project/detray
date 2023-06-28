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

namespace detray::detail {

/// Functors to be used in the @c surface class
template <typename algebra_t>
struct surface_kernels {

    using transform3 = algebra_t;
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

    /// A functor to perform global to local transformation
    struct global_to_local {
        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point3& pos,
            const vector3& dir) const {

            return mask_group[index].to_local_frame(trf3, pos, dir);
        }
    };

    /// A functor to perform local to global transformation
    struct local_to_global {

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline point3 operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3& trf3, const point3& local) const {

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
            // Get the concrete mask instance for this surface
            const auto& m = mask_group[index];
            // Call shape-specific behaviour on the mask
            return m.local_frame().free_to_bound_vector(trf3, free_vec);
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

            const auto& m = mask_group[index];

            return m.local_frame().free_to_bound_jacobian(trf3, free_vec);
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

            const auto& m = mask_group[index];

            return m.local_frame().path_correction(pos, dir, dtds, trf3);
        }
    };
};

}  // namespace detray::detail
