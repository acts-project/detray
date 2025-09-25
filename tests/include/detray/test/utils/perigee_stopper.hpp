/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/constrained_step.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/curvilinear_frame.hpp"

namespace detray {

/// Actor that exits from a navigation stream, if it has reached the perigee
/// @Note Currently only works with step constraints (will be changed in the
/// future)
template <concepts::algebra algebra_t>
struct perigee_stopper : actor {

    using scalar_t = dscalar<algebra_t>;

    struct state {
        /// Radius around the beamline, where to test against the perigee
        /// Outside this radius, the track is considered as not originating
        /// from the IP
        scalar_t m_stopping_radius{10.f * unit<scalar_t>::mm};
        /// Tolerance under which to consider the track at the perigee
        // @TODO Make smaller once overstepping is solved
        scalar_t m_on_perigee_tol{100.f * unit<scalar_t>::um};
        /// Index of the innermost volume for this detector: Convention is 0
        unsigned int m_inner_vol_idx{0u};
    };

    /// Intersects a linear track approximation with the perigee and exits
    /// the navigation, if the perigee is reached.
    ///
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &actor_state,
                                       propagator_state_t &prop_state) const {
        using detector_t = typename propagator_state_t::detector_type;
        using perigee_intersector_t = ray_intersector<line_circular, algebra_t>;

        // Nothing left to do. Propagation will exit successfully on its own
        auto &navigation = prop_state._navigation;
        if (navigation.is_complete()) {
            return;
        }

        // Only check this in the innermost volume and during backward nav.
        if (navigation.volume() != actor_state.m_inner_vol_idx ||
            navigation.direction() != navigation::direction::e_backward) {
            return;
        }

        // Volume that contains the IP
        const tracking_volume inner_vol{navigation.detector(),
                                        actor_state.m_inner_vol_idx};

        // Stop at the perigee (cylindrical detectors only)
        if (inner_vol.id() != volume_id::e_cylinder) {
            return;
        }

        auto &stepping = prop_state._stepping;
        auto &track = stepping();

        // Linear track approximation
        const detail::ray<algebra_t> trk_approx{track.pos(),
                                                -1.f * track.dir()};

        // Check the stopping radius (outside the track will not be stopped)
        assert(actor_state.m_stopping_radius > 0.f);
        constexpr scalar_t max_hz{detail::invalid_value<scalar_t>()};

        const mask<line_circular, algebra_t> perigee_mask{
            actor_state.m_inner_vol_idx, actor_state.m_stopping_radius, max_hz};

        // The perigee is not linked to any detector surface
        constexpr typename detector_t::surface_type inv_sf{};
        const dtransform3D<algebra_t> identity{};
        constexpr scalar_t mask_tolerance{0.f};
        constexpr scalar_t overstep_tolerance{
            -detail::invalid_value<scalar_t>()};

        const auto perigee_intr =
            perigee_intersector_t{}(trk_approx, inv_sf, perigee_mask, identity,
                                    mask_tolerance, overstep_tolerance);

        scalar_t dist_to_cand{std::as_const(navigation).target().path()};
        if (perigee_intr.is_probably_inside() &&
            math::fabs(perigee_intr.path()) < math::fabs(dist_to_cand)) {
            // The track has reached the perigee: "exit success"
            assert(actor_state.m_on_perigee_tol > 0.f);
            if (math::fabs(perigee_intr.path()) <=
                actor_state.m_on_perigee_tol) {
                const curvilinear_frame<algebra_t> cf(track);

                // @TODO: Transport covariance as well
                // assert(!cf.m_bound_vec.is_invalid());
                stepping.bound_params().set_parameter_vector(cf.m_bound_vec);

                prop_state._heartbeat &= navigation.exit();
            } else {
                // @TODO: Use a guided navigator for this in order to catch
                // overstepping correctly
                stepping.template set_constraint<step::constraint::e_actor>(
                    perigee_intr.path());
            }
        }
    }
};

}  // namespace detray
