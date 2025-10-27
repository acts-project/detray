/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/detail/navigation_functions.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"
#include "detray/utils/log.hpp"

namespace detray {

/// @brief The geometry navigation class.
///
/// The navigator is initialized around a detector object, but is itself
/// agnostic to the detectors's object/primitive types.
/// Within a detector volume, the navigatior will perform a local navigation
/// based on the geometry acceleration structure(s) that are provided by the
/// volume. Once the local navigation is resolved, it moves to the next volume
/// by a portal.
/// To this end, it requires a link to the [next] navigation volume in every
/// candidate that is computed by intersection from the detector objects:
/// A module surface must link back to its mother volume, while a portal surface
/// links to the next volume in the direction of the track.
///
/// The navigation state is set up by an init() call and then runs an update()
/// after the track state changed.
///
/// @tparam detector_t the detector to navigate
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t candidate type
template <
    typename detector_t, typename inspector_t = navigation::void_inspector,
    typename intersection_t = intersection2D<typename detector_t::surface_type,
                                             typename detector_t::algebra_type,
                                             !intersection::contains_pos>>
class navigator {

    public:
    using detector_type = detector_t;
    using context_type = detector_type::geometry_context;

    using algebra_type = typename detector_type::algebra_type;
    using scalar_type = dscalar<algebra_type>;
    using point3_type = dpoint3D<algebra_type>;
    using vector3_type = dvector3D<algebra_type>;

    using volume_type = typename detector_type::volume_type;
    using nav_link_type = typename detector_type::surface_type::navigation_link;
    using intersection_type = intersection_t;
    using inspector_type = inspector_t;

    /// @brief Navigation state that contains a the current and the next
    /// candidates
    class state
        : public navigation::base_state<state, detector_type, 2u,
                                        inspector_type, intersection_type> {
        friend class navigator;

        // Allow the filling/updating of candidates
        friend struct intersection_initialize<ray_intersector>;
        friend struct intersection_update<ray_intersector>;

        // Navigation utility functions that need to modify the state
        friend struct navigation::candidate_search;

        template <typename state_t>
        friend constexpr void navigation::update_status(
            state_t &, const navigation::config &);

        template <typename track_t, typename state_t, typename ctx_t>
        friend constexpr void navigation::local_navigation(
            const track_t &, state_t &, const navigation::config &,
            const ctx_t &, const bool);

        template <typename track_t, typename state_t, typename ctx_t>
        friend constexpr void navigation::volume_switch(
            const track_t &, state_t &, const navigation::config &,
            const ctx_t &);

        template <typename track_t, typename state_t, typename ctx_t>
        friend constexpr void navigation::init_loose_cfg(const track_t &,
                                                         state_t &,
                                                         navigation::config,
                                                         const ctx_t &);
        using base_type =
            navigation::base_state<state, detector_type, 2u, inspector_type,
                                   intersection_type>;

        // Result of a geometry object evaluation
        using candidate_t = typename base_type::candidate_t;
        using candidate_cache_t = typename base_type::candidate_cache_t;
        using candidate_itr_t = typename base_type::candidate_itr_t;
        using candidate_const_itr_t = typename base_type::candidate_const_itr_t;
        using dist_t = typename base_type::dist_t;

        public:
        using value_type = candidate_t;

        using view_type = detail::get_view_t<inspector_t>;
        using const_view_type = detail::get_view_t<const inspector_t>;

        /// Use common methods of contructing a nvaigation state
        using base_type::base_type;

        /// @returns the navigator always has only one candidate
        DETRAY_HOST_DEVICE
        constexpr auto n_candidates() const -> dindex { return 1u; }

        /// Update navigation trust level to high trust
        /// @{
        DETRAY_HOST_DEVICE
        constexpr void set_high_trust() {
            this->trust_level(this->trust_level() <
                                      navigation::trust_level::e_high
                                  ? this->trust_level()
                                  : navigation::trust_level::e_high);
        }
        DETRAY_HOST_DEVICE
        constexpr void set_fair_trust() { this->set_high_trust(); }
        /// @}

        private:
        /// Insert a new element @param new_cadidate before position @param pos
        DETRAY_HOST_DEVICE
        constexpr void insert(candidate_const_itr_t /*pos*/,
                              const intersection_type &new_cadidate) {
            // Insert the first candidate
            if (math::fabs(new_cadidate) <
                math::fabs(this->candidates()[1].path)) {
                this->candidates()[1] = new_cadidate;
            }

            assert(this->next_index() == 1u);
            assert(detail::is_invalid_value(new_cadidate.volume_link) ||
                   new_cadidate.volume_link <
                       this->detector().volumes().size());
        }

        /// If the current target is reached, move it to the first position in
        /// the cache and reset the target slot
        DETRAY_HOST_DEVICE
        constexpr void reset_candidate() {
            // Move the current candidate to the first cache position
            this->candidates()[0] = this->candidates()[1];

            // Flag the old candidate as invalid
            this->candidates()[1].set_path(
                std::numeric_limits<scalar_type>::max());
        }
    };

    public:
    /// @brief Helper method to initialize a volume.
    ///
    /// @tparam track_t type of track, needs to provide pos() and dir() methods
    ///
    /// @param track access to the track parameters
    /// @param state the current navigation state
    /// @param cfg the navigation configuration
    /// @param ctx the geometry context
    /// @param resolve_overstepping whether to check for overstepping
    template <typename track_t>
    DETRAY_HOST_DEVICE constexpr void init(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx,
        const bool resolve_overstepping = false) const {

        // Run local navigation in the current volume
        navigation::local_navigation(track, navigation, cfg, ctx,
                                     resolve_overstepping);
    }

    /// @brief Complete update of the navigation flow.
    ///
    /// Restores 'full trust' state to the cadidates cache and checks whether
    /// the track stepped onto a portal and a volume switch is due. If so, or
    /// when the previous update according to the given trust level
    /// failed to restore trust, it performs a complete reinitialization of the
    /// navigation.
    ///
    /// @tparam track_t type of track, needs to provide pos() and dir() methods
    ///
    /// @param track access to the track parameters
    /// @param state the current navigation state
    /// @param cfg the navigation configuration
    /// @param ctx the geometry context
    ///
    /// @returns a heartbeat to indicate if the navigation is still alive
    template <typename track_t>
    DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool update(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx = {},
        const bool /*is_before_actor*/ = true) const {

        assert(navigation.is_alive());
        assert(!track.is_invalid());

        // Update was completely successful (most likely case)
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            DETRAY_VERBOSE_HOST_DEVICE(
                "Nothing left to do (full trust): dist to next %fmm",
                navigation());
            return false;
        }

        // Candidates are re-evaluated based on the current trust level.
        // Should result in 'full trust'
        bool is_init = update_impl(track, navigation, cfg, ctx);

        // Otherwise: if we encountered a portal, perform volume switch
        if (navigation.is_on_portal()) {
            navigation::volume_switch(track, navigation, cfg, ctx);
            is_init = true;

            // Reached end of detector
            if (!navigation.is_alive()) {
                return false;
            }
        }
        // If no trust could be restored during the update, try to rescue the
        // navigation stream by re-initializing with loose tolerances
        if (navigation.trust_level() != navigation::trust_level::e_full ||
            navigation.cache_exhausted()) {

            is_init = true;
            navigation::init_loose_cfg(track, navigation, cfg, ctx);

            navigation.run_inspector(cfg, track.pos(), track.dir(),
                                     "Re-init: ");
        }

        DETRAY_VERBOSE_HOST_DEVICE("Update complete: dist to next %fmm",
                                   navigation());

        return is_init;
    }

    private:
    /// Helper method to update the candidate (surface intersections)
    /// based on an externally provided trust level. Will (re-)initialize the
    /// navigation if there is no trust.
    ///
    /// @tparam track_t type of track, needs to provide pos() and dir() methods
    ///
    /// @param track access to the track parameters
    /// @param state the current navigation state
    /// @param cfg the navigation configuration
    /// @param ctx the geometry context
    template <typename track_t>
    DETRAY_HOST_DEVICE constexpr bool update_impl(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx) const {

        const auto &det = navigation.detector();
        constexpr bool is_init{true};

        assert(navigation.trust_level() != navigation::trust_level::e_full);

        // Update only the current candidate and the corresponding
        // - do this only when the navigation state is still coherent
        if (navigation.trust_level() == navigation::trust_level::e_high) {

            DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - high trust");

            // Update next candidate: If not reachable, 'high trust' is broken
            if (!navigation::update_candidate(
                    navigation.direction(), navigation.target(), track, det,
                    cfg.intersection, navigation.external_tol(), ctx)) {
                navigation.status(navigation::status::e_unknown);
                // This will run into the fair trust case below.
                navigation.set_fair_trust();
            } else {
                // Update navigation flow on the new candidate information
                navigation::update_status(navigation, cfg);

                navigation.run_inspector(cfg, track.pos(), track.dir(),
                                         "Update complete: high trust: ");

                // The work is done if: the track has not reached a surface yet
                // or trust is gone (portal was reached or the cache is broken).
                if (navigation.status() ==
                        navigation::status::e_towards_object ||
                    navigation.is_on_portal()) {
                    return !is_init;
                }
                // Else (if full trust): Track is on non-portal surface and
                // cache is not exhausted. Move the current surface back in the
                // cache and re-initialize the volume to find the next target
                if (navigation.trust_level() ==
                    navigation::trust_level::e_full) {
                    navigation.reset_candidate();
                }

                // Find the next/different candidate
                navigation.set_no_trust();
            }
        }
        // If [next] target is not reachable or actor flagged 'no trust',
        // re-initialize the volume.
        assert(navigation.trust_level() == navigation::trust_level::e_no_trust);
        DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - no trust");

        constexpr bool resolve_overstepping{true};
        navigation::local_navigation(track, navigation, cfg, ctx,
                                     resolve_overstepping);
        return is_init;
    }
};

}  // namespace detray
