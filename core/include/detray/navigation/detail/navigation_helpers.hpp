/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/detail/navigation_base_state.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/tracks/ray.hpp"

namespace detray::navigation {

/// A functor that fills the navigation candidates cache by intersecting
/// a surface in the track position neighbourhood with the tangential to
/// the track direction
struct candidate_search {

    /// Visit the volume acceleration data structure: Returns a range of
    /// surfaces each of which is passed to this functor
    ///
    /// @param sf_descr descriptor of one of the surfaces in the track
    ///                 neighbourhood
    /// @param det access to the detector
    /// @param ctx the geometry context
    /// @param track the track parameters
    /// @param nav_state state of navigation stream of the track
    /// @param mask_tol min. and max. mask tolerance
    /// @param mask_tol_scalor scale factor with track distance in range of
    ///                        @c mask_tol
    /// @param overstep_tol how far behind the track pos to look for
    /// candidates
    template <typename track_t, typename detector_t,
              typename navigation_state_t>
    DETRAY_HOST_DEVICE constexpr void operator()(
        const typename detector_t::surface_type &sf_descr,
        const detector_t &det, const typename detector_t::geometry_context &ctx,
        const track_t &track, navigation_state_t &nav_state,
        const darray<typename detector_t::scalar_type, 2> mask_tol,
        const typename detector_t::scalar_type mask_tol_scalor,
        const typename detector_t::scalar_type overstep_tol) const {

        using algebra_t = typename detector_t::algebra_type;
        using scalar_t = dscalar<algebra_t>;

        const auto sf = tracking_surface{det, sf_descr};

        // Tangential to the track direction
        detail::ray<algebra_t> tangential{
            track.pos(),
            static_cast<scalar_t>(nav_state.direction()) * track.dir()};

        // Perform intersection and add result to the navigation cache via
        // @c nav_state.insert()
        sf.template visit_mask<intersection_initialize<ray_intersector>>(
            nav_state, tangential, sf_descr, det.transform_store(), ctx,
            sf.is_portal() ? darray<scalar_t, 2>{0.f, 0.f} : mask_tol,
            mask_tol_scalor, overstep_tol);
    }
};

/// @returns true if the candidate lies on a surface
template <typename candidate_t>
DETRAY_HOST_DEVICE constexpr bool has_reached_surface(
    const candidate_t &candidate, const navigation::config &cfg) {
    return (math::fabs(candidate.path) < cfg.path_tolerance);
}

/// @brief Helper method that updates the intersection of a single candidate
/// and checks reachability
///
/// @tparam candidate_t type of navigation candidate (intersection result)
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam detector_t type of the detector
/// @tparam context_t type of geometry context
///
/// @param nav_dir the navigation direction
/// @param candidate the candidate intersection to be updated
/// @param track access to the track parameters
/// @param det access to the detector (geometry)
/// @param cfg the navigation configuration
/// @param ctx the geometry context
///
/// @returns @c true if the track can reach this candidate.
template <typename candidate_t, typename track_t, typename detector_t>
DETRAY_HOST_DEVICE inline bool update_candidate(
    const navigation::direction nav_dir, candidate_t &candidate,
    const track_t &track, const detector_t &det, const navigation::config &cfg,
    const typename detector_t::geometry_context &ctx) {

    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    // Invalid intersection result cannot be updated
    if (candidate.sf_desc.barcode().is_invalid()) {
        return false;
    }

    const auto sf = tracking_surface{det, candidate.sf_desc};

    // Tangential to the track direction
    auto tangential{detail::ray<algebra_t>(
        track.pos(), static_cast<scalar_t>(nav_dir) * track.dir())};

    // Minimum and maximum tolerance for mask check
    auto mask_tol{sf.is_portal() ? darray<scalar_t, 2>{0.f, 0.f}
                                 : darray<scalar_t, 2>{cfg.min_mask_tolerance,
                                                       cfg.max_mask_tolerance}};

    // Perform intersection and check whether this candidate is reachable by
    // the track
    return sf.template visit_mask<intersection_update<ray_intersector>>(
        std::move(tangential), candidate, det.transform_store(), ctx, mask_tol,
        static_cast<scalar_t>(cfg.mask_tolerance_scalor),
        static_cast<scalar_t>(cfg.overstep_tolerance));
}

/// @brief Helper method that re-establishes the navigation status after an
/// update.
///
/// It checks wether the track has reached a surface or is still moving
/// towards the next surface candidate. If no new next candidate can be
//  found, it flags 'no trust' in order to trigger a volume initialization.
///
/// @param navigation the current navigation state
/// @param cfg the navigation configuration
template <typename navigation_state_t>
DETRAY_HOST_DEVICE inline void update_status(navigation_state_t &navigation,
                                             const navigation::config &cfg) {

    // Check whether the track reached the current candidate. Might be a
    // portal, in which case the navigation needs to be re-initialized
    if (!navigation.is_exhausted() &&
        navigation::has_reached_surface(navigation.target(), cfg)) {
        // Set the next object that we want to reach (this function is only
        // called once the cache has been updated to a full trust state).
        // Might lead to exhausted cache.
        navigation.advance();
        navigation.set_status(navigation.get_current().sf_desc.is_portal()
                                  ? navigation::status::e_on_portal
                                  : navigation::status::e_on_object);
    } else {
        // Otherwise the track is moving towards a surface
        navigation.set_status(navigation::status::e_towards_object);
    }
    // Exhaustion happens when after an update no next candidate in the
    // cache is reachable anymore -> triggers init of [new] volume
    // In backwards navigation or with strongly bent tracks, the cache may
    // not be exhausted when trying to exit the volume (the ray is seeing
    // the opposite side of the volume)
    navigation.set_trust_level(navigation.is_exhausted() ||
                                       navigation.on_portal()
                                   ? navigation::trust_level::e_no_trust
                                   : navigation::trust_level::e_full);
}

/// @brief Helper method to initialize a navigation state in a given volume.
///
/// Calls the volumes acceleration structure, then tests the surfaces for
/// intersection and keeps the clostest one(s) ("local navigation" in the
/// volume). The closest candidate is set as 'next candidate' or 'target'.
///
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam navigation_state_t the state type of the navigation stream
/// @tparam context_t the type of geometry context
///
/// @param track access to the track parameters
/// @param navigation the current navigation state
/// @param cfg the navigation configuration
/// @param ctx the geometry context
template <typename track_t, typename navigation_state_t, typename context_t>
DETRAY_HOST_DEVICE inline void local_navigation(const track_t &track,
                                                navigation_state_t &navigation,
                                                const navigation::config &cfg,
                                                const context_t &ctx) {
    using algebra_t = typename track_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    const auto &det = navigation.detector();
    const auto volume = tracking_volume{det, navigation.volume()};

    // Clean up state
    navigation.clear_cache();
    navigation.set_heartbeat(true);

    // Search for neighboring surfaces and fill candidates into cache
    volume.template visit_neighborhood<candidate_search>(
        track, cfg, ctx, det, ctx, track, navigation,
        darray<scalar_t, 2u>{cfg.min_mask_tolerance, cfg.max_mask_tolerance},
        static_cast<scalar_t>(cfg.mask_tolerance_scalor),
        static_cast<scalar_t>(cfg.overstep_tolerance));

    // Determine overall status of the navigation after updating the cache
    update_status(navigation, cfg);

    // If init was not successful, the navigation setup is broken
    navigation.set_heartbeat(navigation.trust_level() ==
                             navigation::trust_level::e_full);
}

/// @brief Perform a detector volume switch.
///
/// Once a portal is reached, this function will update the navigation
/// stream to continue in the new volume. If it is a valid detector volume,
/// the navigation is re-initialized by performing local navigation in the
/// volume. If the end of the detector geometry was reached, the navigation
/// exits.
///
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam navigation_state_t the state type of the navigation stream
/// @tparam context_t the type of geometry context
///
/// @param track access to the track parameters
/// @param navigation the current navigation state
/// @param cfg the navigation configuration
/// @param ctx the geometry context
///
/// @returns a boolean that indicates wheather the volume was
/// (re-)initialized
template <typename track_t, typename navigation_state_t, typename context_t>
DETRAY_HOST_DEVICE inline bool volume_switch(const track_t &track,
                                             navigation_state_t &navigation,
                                             const navigation::config &cfg,
                                             const context_t &ctx) {

    // Set volume index to the next volume provided by the portal
    navigation.set_volume(navigation.get_current().volume_link);

    // Navigation reached the end of the detector world
    if (detail::is_invalid_value(navigation.volume())) {
        navigation.exit();
        navigation.run_inspector(cfg, track.pos(), track.dir(),
                                 "Reached end of detector: ");
        return false;
    }

    // Either end of world or valid volume index
    assert(navigation.volume() < navigation.detector().volumes().size());

    // Init the new volume
    local_navigation(track, navigation, cfg, ctx);

    // Fresh initialization: reset trust and hearbeat even though we are
    // on inner portal
    navigation.set_trust_level(navigation::trust_level::e_full);
    navigation.set_heartbeat(!navigation.is_exhausted());

    if (!navigation.is_alive()) {
        navigation.abort();
    }

    return true;
}

/// @brief Initilaize the volume with loose configuration.
///
/// If trust cannot be established and/or no surfaces can be found in the
/// current volume anymore, try to save the navigation stream by looking
/// for candidates further behind the track
///
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam navigation_state_t the state type of the navigation stream
/// @tparam context_t the type of geometry context
///
/// @param track access to the track parameters
/// @param navigation the current navigation state
/// @param cfg the navigation configuration (copy on function stack)
/// @param ctx the geometry context
///
/// @returns a boolean that indicates wheather the volume was
/// (re-)initialized
template <typename track_t, typename navigation_state_t, typename context_t>
DETRAY_HOST_DEVICE inline bool init_loose_cfg(const track_t &track,
                                              navigation_state_t &navigation,
                                              navigation::config loose_cfg,
                                              const context_t &ctx) {

    // Use the max mask tolerance in case a track leaves the volume
    // when a sf is 'sticking' out of the portals due to the tol
    const auto new_overstep_tol{
        math::min(100.f * loose_cfg.overstep_tolerance,
                  -10.f * loose_cfg.max_mask_tolerance)};
    loose_cfg.overstep_tolerance = new_overstep_tol;

    local_navigation(track, navigation, loose_cfg, ctx);

    // Unrecoverable
    if (navigation.trust_level() != navigation::trust_level::e_full) {
        navigation.abort();
    }

    return true;
}

}  // namespace detray::navigation
