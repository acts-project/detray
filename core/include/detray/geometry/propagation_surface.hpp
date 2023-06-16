
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detail/propagation_kernels.hpp"
#include "detray/geometry/surface.hpp"

namespace detray::propagation {

/// @brief Facade for a detray detector surface with propagation functionality.
///
/// Extends a geometryic detector surface with geometry dependent propagation
/// functionality like bound-to-free Jacobian generation.
template <typename detector_t>  // @TODO: This needs a concept
class surface : public detray::surface<detector_t> {

    using base_type = detray::surface<detector_t>;

    /// Surface descriptor type
    using descr_t = typename detector_t::surface_type;

    public:
    using algebra = typename base_type::transform3;
    using transform3 = algebra;
    using point3 = typename base_type::point3;
    using vector3 = typename base_type::point3;
    using context = typename base_type::context;

    /// Vector type for track parameters in global coordinates
    using free_vector_type =
        typename free_track_parameters<algebra>::vector_type;
    /// Vector type for track parameters in local (bound) coordinates
    using bound_vector_type =
        typename bound_track_parameters<algebra>::vector_type;

    // Import base class constructors
    using base_type::base_type;

    /// @returns the track parametrization projected onto the surface (bound)
    DETRAY_HOST_DEVICE
    constexpr auto free_to_bound_vector(
        const context &ctx, const free_vector_type &free_vec) const {
        return this->template visit_mask<detail::free_to_bound_vector<algebra>>(
            this->transform(ctx), free_vec);
    }

    /// @returns the global track parametrization from a bound representation
    DETRAY_HOST_DEVICE
    constexpr auto bound_to_free_vector(
        const context &ctx, const bound_vector_type &bound_vec) const {
        return this->template visit_mask<detail::bound_to_free_vector<algebra>>(
            this->transform(ctx), bound_vec);
    }

    /// @returns the jacobian to go from a free to a bound track parametrization
    DETRAY_HOST_DEVICE
    constexpr auto free_to_bound_jacobian(
        const context &ctx, const free_vector_type &free_vec) const {
        return this
            ->template visit_mask<detail::free_to_bound_jacobian<algebra>>(
                this->transform(ctx), free_vec);
    }

    /// @returns the jacobian to go from a bound to a free track parametrization
    DETRAY_HOST_DEVICE
    constexpr auto bound_to_free_jacobian(
        const context &ctx, const bound_vector_type &bound_vec) const {
        return this
            ->template visit_mask<detail::bound_to_free_jacobian<algebra>>(
                this->transform(ctx), bound_vec);
    }

    /// @returns the path correction term
    DETRAY_HOST_DEVICE
    constexpr auto path_correction(const context &ctx, const vector3 &pos,
                                   const vector3 &dir,
                                   const vector3 &dtds) const {
        return this->template visit_mask<detail::path_correction<algebra>>(
            this->transform(ctx), pos, dir, dtds);
    }
};

template <typename detector_t, typename descr_t>
DETRAY_HOST_DEVICE surface(const detector_t &, const descr_t &)
    ->surface<detector_t>;

template <typename detector_t>
DETRAY_HOST_DEVICE surface(const detector_t &, const geometry::barcode)
    ->surface<detector_t>;

}  // namespace detray::propagation
