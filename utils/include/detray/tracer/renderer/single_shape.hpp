/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/propagator/base_actor.hpp"

// System include(s)
#include <iostream>
#include <limits>
#include <memory>

namespace detray {

/// Calculates the color of a pixel. Starting point of the shader pipeline
template <typename T, template <typename> class algebra_t, typename mask_t,
          typename material_t = material_slab<T>>
struct single_shape : detray::actor {

    using link_t = dsimd<algebra_t, std::uint_least16_t>;
    using surface_t = surface_descriptor<dtyped_index<>, dtyped_index<>, link_t,
                                         link_t, link_t>;
    using intersection_t = intersection2D<surface_t, T, algebra_t>;

    struct global_state {
        using scalar_t = dscalar<algebra_t<T>>;
        using point3D = dpoint3D<algebra_t<T>>;
        using vector3D = dvector3D<algebra_t<T>>;
        using transform3D = dtransform3D<algebra_t<T>>;

        /// Construct from surface data:
        DETRAY_HOST_DEVICE
        global_state(const transform3D &trf, const mask_t &mask,
                     const material_t &mat)
            : m_trf{trf}, m_mask{mask}, m_material{mat} {}

        /// Threadsafe interface
        /// @{
        const transform3D &transform() const { return m_trf; }
        const mask_t &mask() const { return m_mask; }
        const material_t &material() const { return m_material; }
        /// @}

        /// The surfaces data
        transform3D m_trf;
        mask_t m_mask;
        material_t m_material;
    };

    struct state {

        DETRAY_HOST_DEVICE
        state(const T min = 0.f, const T max = std::numeric_limits<T>::max())
            : m_interval{min, max} {}

        const material_t &material() const { return *m_material; }

        /// From potentialliy multiple intersected surfaces in the intersection,
        /// get the index of the closest one
        std::size_t closest_solution() const {
            return m_intersections[0].status.firstOne();
        }

        std::array<T, 2> m_interval;
        /// Resulting intersection
        std::array<intersection_t, 2> m_intersections{};
        /// Pointer to the material of the surface
        const material_t *m_material{nullptr};
        /// Flag to the obseving colorizer/shaders that the surface was hit
        bool m_is_inside = false;
    };

    /// Intersect the ray with the mask. The closest intersection is in front of
    /// the @c m_intersections container
    template <typename scene_handle_t>
    DETRAY_HOST_DEVICE void operator()(state &loc_st,
                                       const scene_handle_t &sc) const {
        // In this special case, the geometry will be this actor's global_state
        const global_state &geo = sc.geometry();

        // Perform the intersection
        loc_st.m_is_inside = place_in_collection(
            geo.mask().template intersector<intersection_t>()(
                sc.ray(), surface_t{}, geo.mask(), geo.transform()),
            loc_st.m_intersections);
        if (loc_st.m_is_inside) {
            loc_st.m_material = std::addressof(geo.material());
        }
    }

    private:
    /// Places the single solution of a ray-surface intersection @param sfi
    /// in the given container @param intersections, if the surfaces was hit.
    /// @todo: Add sorting algorithm
    ///
    /// @returns @c true if the intersection was is valid.
    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        typename is_container_t::value_type &&sfi,
        is_container_t &intersections) const {
        if (sfi.is_inside()) {
            intersections[0] = sfi;
            return true;
        } else {
            return false;
        }
    }

    /// Places all of those solutions of a ray-surface intersection @param sfi
    /// in the given container @param intersections, that hit the surface
    /// @todo: Add sorting algorithm
    ///
    /// @returns @c true if at least one valid intersection solution was found.
    template <typename is_container_t>
    DETRAY_HOST_DEVICE bool place_in_collection(
        std::array<typename is_container_t::value_type, 2> &&solutions,
        is_container_t &intersections) const {
        bool is_valid = false;
        std::size_t n_sol{0u};
        for (std::size_t i = 0u; i < 2u; ++i) {
            if (solutions[i].is_inside()) {
                intersections[n_sol++] = solutions[i];
                is_valid = true;
            }
        }
        return is_valid;
    }
};

}  // namespace detray