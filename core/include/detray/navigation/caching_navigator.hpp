/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/navigation_functions.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_kernel.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"
#include "detray/utils/log.hpp"

namespace detray {

namespace navigation {

static constexpr std::size_t default_cache_size{8u};

}  // namespace navigation

/// @brief Navigation class which caches a 'road' through the detector
///
/// This navigator applies a trust level based update of its candidate
/// (intersection) cache, which is kept in the naviagtor's state. The trust
/// level, and with it the appropriate update policy, must be set by an actor,
/// otherwise no update will be performed.
///
/// @tparam detector_t the detector to navigate
/// @tparam k_cache_capacity the capacity of the candidate cache
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t candidate type
template <typename detector_t,
          std::size_t k_cache_capacity = navigation::default_cache_size,
          typename inspector_t = navigation::void_inspector,
          typename intersection_t = intersection2D<
              typename detector_t::surface_type,
              typename detector_t::algebra_type, !intersection::contains_pos>>
class caching_navigator {

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

    /// @brief Navigation state that contains a cache of candidates
    ///
    /// Once a volume is reached, the cache is filled by building a 'road' of
    /// surfaces that will likely be encountered by the track in that volume.
    /// The cache keeps a range of reachable candidates that lie between the
    /// next and last index in the cache. The cache is sorted by distance
    /// to the track position.
    class state
        : public navigation::base_state<state, detector_type, k_cache_capacity,
                                        inspector_type, intersection_type> {

        friend class caching_navigator;

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
            navigation::base_state<state, detector_type, k_cache_capacity,
                                   inspector_type, intersection_type>;

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

        /// Update navigation trust level to high trust
        DETRAY_HOST_DEVICE
        constexpr void set_high_trust() {
            this->trust_level(this->trust_level() <
                                      navigation::trust_level::e_high
                                  ? this->trust_level()
                                  : navigation::trust_level::e_high);
        }

        /// Update navigation trust level to fair trust
        DETRAY_HOST_DEVICE
        constexpr void set_fair_trust() {
            this->trust_level(this->trust_level() <
                                      navigation::trust_level::e_fair
                                  ? this->trust_level()
                                  : navigation::trust_level::e_fair);
        }

        private:
        /// Insert a new element @param new_candidate before position @param pos
        DETRAY_HOST_DEVICE
        constexpr void insert(candidate_const_itr_t pos,
                              const intersection_type &new_candidate) {

            // Candidate is too far away to be placed in cache
            if (pos == this->candidates().end()) {
                return;
            }

            assert(detail::is_invalid_value(new_candidate.volume_link) ||
                   new_candidate.volume_link <
                       this->detector().volumes().size());

            // Insert the first candidate
            if (this->n_candidates() == 0) {
                this->candidates()[0] = new_candidate;
                this->last_index(this->last_index() + 1);
                assert(this->next_index() <= this->last_index() + 1);
                assert(static_cast<std::size_t>(this->last_index()) <
                       k_cache_capacity);
                return;
            }

            // Position where to insert the new candidate
            auto idx{static_cast<dist_t>(
                detray::ranges::distance(this->candidates().cbegin(), pos))};
            assert(idx >= 0);

            // Do not add the same surface (intersection) multiple times
            const auto is_clash_at_pos = [this,
                                          &new_candidate](std::size_t index) {
                return (this->candidates()[index].sf_desc.barcode() ==
                        new_candidate.sf_desc.barcode()) &&
                       (math::fabs(this->candidates()[index].path() -
                                   new_candidate.path()) <=
                        1.f * unit<scalar_type>::um);
            };

            const auto idxu{static_cast<std::size_t>(idx)};
            if (is_clash_at_pos(idxu) ||
                ((idxu > 0u) && is_clash_at_pos(idxu - 1u))) {
                return;
            }

            // Shift all following candidates and evict the last element,
            // if the cache is already full
            constexpr auto shift_max{static_cast<dist_t>(k_cache_capacity - 2)};
            const dist_t shift_begin{math::min(this->last_index(), shift_max)};

            for (dist_t i = shift_begin; i >= idx; --i) {
                const auto j{static_cast<std::size_t>(i)};
                this->candidates()[j + 1u] = this->candidates()[j];
            }

            // Now insert the new candidate and update candidate range
            this->candidates()[static_cast<std::size_t>(idx)] = new_candidate;
            this->last_index(
                math::min(static_cast<dist_t>(this->last_index() + 1),
                          static_cast<dist_t>(k_cache_capacity - 1)));

            assert(this->next_index() <= this->last_index() + 1);
            assert(static_cast<std::size_t>(this->last_index()) <
                   k_cache_capacity);
        }

        /// Clear the state
        DETRAY_HOST_DEVICE constexpr void clear_cache() {
            base_type::clear_cache();
            this->next_index(0);
            this->last_index(-1);
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

        // Candidates are re-evaluated based on the current trust level.
        // Should result in 'full trust'
        bool is_init = update_impl(track, navigation, cfg, ctx);

        // Update was completely successful (most likely case)
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            DETRAY_VERBOSE_HOST_DEVICE("Update complete: dist to next %fmm",
                                       navigation());
            return is_init;
        }
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
    /// Helper method to update the candidates (surface intersections)
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

        // Current candidates are up to date, nothing left to do
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - full trust");
            return !is_init;
        }
        // Update only the current candidate and the corresponding next target
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
                // cache is not exhausted. Ready the next target
                if (navigation.trust_level() ==
                        navigation::trust_level::e_full &&
                    navigation::update_candidate(
                        navigation.direction(), navigation.target(), track, det,
                        cfg.intersection, navigation.external_tol(), ctx)) {
                    return !is_init;
                }

                // If next candidate is not reachable, don't 'return', but
                // escalate the trust level.
                // This will run into the fair trust case below or the no trust
                // case if the cache is broken
                navigation.set_fair_trust();
            }
        }
        // Re-evaluate all currently available candidates and sort again
        // - do this when your navigation state is stale, but not invalid
        if (navigation.trust_level() == navigation::trust_level::e_fair &&
            !navigation.cache_exhausted()) {

            DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - fair trust");

            for (auto &candidate : navigation) {
                // Disregard this candidate if it is not reachable
                if (!navigation::update_candidate(
                        navigation.direction(), candidate, track, det,
                        cfg.intersection, navigation.external_tol(), ctx)) {
                    // Forcefully set dist to numeric max for sorting
                    candidate.set_path(std::numeric_limits<scalar_type>::max());
                }
            }
            detray::sequential_sort(navigation.begin(), navigation.end());
            // Take the nearest (sorted) candidate first
            navigation.set_next(navigation.begin());
            // Ignore unreachable elements (needed to determine exhaustion)
            navigation.set_last(find_invalid(navigation.candidates()));
            // Update navigation flow on the new candidate information
            navigation::update_status(navigation, cfg);

            navigation.run_inspector(cfg, track.pos(), track.dir(),
                                     "Update complete: fair trust: ");

            // If there are no reachable candidates in the cache after
            // re-evaluation, re-initialize the volume
            if (navigation.cache_exhausted()) {
                navigation.set_no_trust();
            }
        }
        // Re-initialize the volume (actor flagged 'no trust' or previous trust
        // level update failed)
        if (navigation.trust_level() == navigation::trust_level::e_no_trust) {

            DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - no trust");

            constexpr bool resolve_overstepping{true};
            navigation::local_navigation(track, navigation, cfg, ctx,
                                         resolve_overstepping);
            return is_init;
        }

        return !is_init;
    }

    /// Helper to evict all unreachable/invalid candidates from the cache:
    /// Finds the first unreachable candidate (has been invalidated during
    /// update) in a sorted (!) cache.
    ///
    /// @param candidates the cache of candidates to be cleaned
    ///
    /// @returns iterator to the last reachable candidate
    DETRAY_HOST_DEVICE constexpr auto find_invalid(
        const typename state::candidate_cache_t &candidates) const {

        // Depends on previous invalidation of unreachable candidates!
        auto not_reachable = [](const intersection_type &candidate) {
            return candidate.path() == std::numeric_limits<scalar_type>::max();
        };

        return detray::find_if(candidates.begin(), candidates.end(),
                               not_reachable);
    }
};

}  // namespace detray
