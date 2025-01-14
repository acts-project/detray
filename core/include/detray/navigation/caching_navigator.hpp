/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/navigation/detail/navigation_base_state.hpp"
#include "detray/navigation/detail/navigation_helpers.hpp"
#include "detray/navigation/detail/print_state.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <ostream>

namespace detray {

namespace navigation {

/// Default size of the internal navigation state candaidate cache
static constexpr std::size_t default_cache_size{10u};

}  // namespace navigation

/// @brief A navigation class that performs "road building" per volume.
///
/// This navigator applies a trust level based update of its candidate
/// (intersection) cache, which is kept in the naviagtor's state. The trust
/// level, and with it the appropriate update policy, must be set by an actor,
/// otherwise no update will be performed.
/// The cache contains a "road" of surfaces within the current volume, which is
/// comprised of the surfaces that are most likely to be encountered by the
/// track as it is propagated through the volume. After a volume switch, the
/// road is re-built. During volume traversal, the candidates are replayed
/// from the cache.
///
/// @tparam detector_t the detector to navigate
/// @tparam k_cache_capacity the capacity of the candidate cache
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t candidate type
template <typename detector_t,
          std::size_t k_cache_capacity = navigation::default_cache_size,
          typename inspector_t = navigation::void_inspector,
          typename intersection_t =
              intersection2D<typename detector_t::surface_type,
                             typename detector_t::algebra_type, false>>
class caching_navigator {

    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    using context_t = detector_t::geometry_context;

    public:
    using detector_type = detector_t;
    using intersection_type = intersection_t;
    using inspector_type = inspector_t;

    /// @brief Navigation state that contains a cache of candidates
    ///
    /// The cache keeps a range of reachable candidates that lie between the
    /// @c m_next and @c m_last index in the cache. The cache is sorted
    /// by distance from the track position. It is updated by the navigator
    /// according to different "trust levels" which model the degree to which
    /// the data in the cache has become stale after updates to the track params
    class state final
        : public navigation::base_state<detector_t, k_cache_capacity,
                                        inspector_t, intersection_t> {

        using base_t = navigation::base_state<detector_t, k_cache_capacity,
                                              inspector_t, intersection_t>;

        using point3_t = dpoint3D<algebra_t>;
        using vector3_t = dvector3D<algebra_t>;

        using candidate_itr_t = typename base_t::candidate_itr_t;
        using candidate_const_itr_t = typename base_t::candidate_const_itr_t;
        using dist_t = typename base_t::dist_t;

        public:
        // The navigator can modify the cache
        friend class caching_navigator;

        // Navigation utility functions that need to modify the state
        friend struct candidate_search;
        template <typename state_t>
        friend void navigation::update_status(state_t &,
                                              const navigation::config &);
        template <typename track_t, typename state_t, typename ctx_t>
        friend void navigation::local_navigation(const track_t &, state_t &,
                                                 const navigation::config &,
                                                 const ctx_t &);
        template <typename track_t, typename state_t, typename ctx_t>
        friend bool navigation::volume_switch(const track_t &, state_t &,
                                              const navigation::config &,
                                              const ctx_t &);
        template <typename track_t, typename state_t, typename ctx_t>
        friend bool navigation::init_loose_cfg(const track_t &, state_t &,
                                               navigation::config,
                                               const ctx_t &);

        // Allow the filling/updating of candidates
        friend struct intersection_initialize<ray_intersector>;
        friend struct intersection_update<ray_intersector>;

        // Need at least two slots for the caching to work
        static_assert(
            k_cache_capacity >= 2u,
            "Navigation cache needs to have a capacity larger than 1");

        /// Use common methods of contructing a nvaigation state
        using base_t::base_t;

        private:
        /// Set the next surface that we want to reach (update target)
        DETRAY_HOST_DEVICE
        constexpr void advance() {
            this->set_next_index(this->next_index() + 1);
        }

        /// Set the next surface that we want to reach (update target)
        DETRAY_HOST_DEVICE
        constexpr void set_next(candidate_itr_t new_next) {
            const auto new_idx{static_cast<dist_t>(detray::ranges::distance(
                this->get_candidates().begin(), new_next))};
            this->set_next_index(new_idx);
        }

        /// Updates the position of the last valid candidate
        DETRAY_HOST_DEVICE
        constexpr void set_last(candidate_const_itr_t new_last) {
            const auto new_idx{static_cast<dist_t>(
                detray::ranges::distance(static_cast<candidate_const_itr_t>(
                                             this->get_candidates().begin()),
                                         new_last) -
                1)};
            this->set_last_index(new_idx);
        }

        /// Insert a new element @param new_cadidate before position @param pos
        /// This is the cache filling scheme that is called during
        /// 'local_navigation'
        DETRAY_HOST_DEVICE
        constexpr void insert(candidate_const_itr_t pos,
                              const intersection_type &new_cadidate) {

            // Candidate is too far away to be placed in cache
            if (pos == this->get_candidates().end()) {
                return;
            }

            assert(detail::is_invalid_value(new_cadidate.volume_link) ||
                   new_cadidate.volume_link <
                       this->detector().volumes().size());

            // Insert the first candidate
            if (this->is_exhausted()) {
                this->get_candidates()[0] = new_cadidate;
                this->set_last_index(this->last_index() + 1);
                assert(this->next_index() <= this->last_index() + 1);
                assert(static_cast<std::size_t>(this->last_index()) <
                       k_cache_capacity);
                return;
            }

            // Position where to insert the new candidate
            auto idx{static_cast<dist_t>(detray::ranges::distance(
                std::as_const(*this).get_begin(), pos))};
            assert(idx >= 0);

            // Shift all following candidates and evict the last element,
            // if the cache is already full
            constexpr auto shift_max{static_cast<dist_t>(k_cache_capacity - 2)};
            const dist_t shift_begin{math::min(this->last_index(), shift_max)};

            for (dist_t i = shift_begin; i >= idx; --i) {
                const auto j{static_cast<std::size_t>(i)};
                this->get_candidates()[j + 1u] = this->get_candidates()[j];
            }

            // Now insert the new candidate and update candidate range
            this->get_candidates()[static_cast<std::size_t>(idx)] =
                new_cadidate;
            this->set_last_index(
                math::min(static_cast<dist_t>(this->last_index() + 1),
                          static_cast<dist_t>(k_cache_capacity - 1)));

            assert(this->next_index() <= this->last_index() + 1);
            assert(static_cast<std::size_t>(this->last_index()) <
                   k_cache_capacity);
        }

        /// Call the navigation inspector
        DETRAY_HOST_DEVICE
        inline void run_inspector(
            [[maybe_unused]] const navigation::config &cfg,
            [[maybe_unused]] const point3_t &track_pos,
            [[maybe_unused]] const vector3_t &track_dir,
            [[maybe_unused]] const char *message) {
            if constexpr (!std::is_same_v<inspector_t,
                                          navigation::void_inspector>) {
                this->get_inspector()(*this, cfg, track_pos, track_dir,
                                      message);
            }
        }

        private:
        /// Transform to a string for debugging output
        DETRAY_HOST
        friend std::ostream &operator<<(std::ostream &out_stream,
                                        const state &state) {
            out_stream << navigation::print_state(state);

            return out_stream;
        }
    };

    public:
    /// @brief Initialize the navigation stream: Start of navigation
    ///
    /// @tparam track_t type of track, needs to provide pos() and dir() methods
    ///
    /// @param track access to the track parameters
    /// @param navigation the current navigation state
    /// @param cfg the navigation configuration
    /// @param ctx the geometry context
    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(const track_t &track,
                                        caching_navigator::state &navigation,
                                        const navigation::config &cfg,
                                        const context_t &ctx) const {
        // Run local navigation in the current volume
        navigation::local_navigation(track, navigation, cfg, ctx);

        navigation.run_inspector(cfg, track.pos(), track.dir(),
                                 "Init complete: ");
    }

    /// @brief Complete update of the navigation stream.
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
    /// @param navigation the current navigation state
    /// @param cfg the navigation configuration
    /// @param ctx the geometry context
    ///
    /// @returns a heartbeat to indicate if the navigation is still alive
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update(const track_t &track,
                                          caching_navigator::state &navigation,
                                          const navigation::config &cfg,
                                          const context_t &ctx = {}) const {
        // Candidates are re-evaluated based on the current trust level.
        // Should result in 'full trust'
        bool is_init = update_impl(track, navigation, cfg, ctx);

        // Update was completely successful (most likely case)
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            return is_init;
        }
        // Otherwise: did we run into a portal?
        else if (navigation.on_portal()) {
            navigation::volume_switch(track, navigation, cfg, ctx);

            navigation.run_inspector(cfg, track.pos(), track.dir(),
                                     "Volume switch: ");
        }
        // If no trust could be restored for the current state, (local)
        // navigation might be exhausted: re-initialize volume
        else {
            init(track, navigation, cfg, ctx);
            is_init = true;

            // Sanity check: Should never be the case after complete update call
            if (navigation.trust_level() != navigation::trust_level::e_full) {
                navigation::init_loose_cfg(track, navigation, cfg, ctx);
            }

            navigation.run_inspector(cfg, track.pos(), track.dir(),
                                     "Re-init: ");
        }

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
    /// @param navigation the current navigation state
    /// @param cfg the navigation configuration
    /// @param ctx the geometry context
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update_impl(
        const track_t &track, caching_navigator::state &navigation,
        const navigation::config &cfg, const context_t &ctx) const {

        const auto &det = navigation.detector();

        // Current candidates are up to date, nothing left to do
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            return false;
        }

        // Update only the current candidate and the corresponding next target
        // - do this only when the navigation state is still coherent
        if (navigation.trust_level() == navigation::trust_level::e_high) {
            // Update next candidate: If not reachable, 'high trust' is broken
            if (!navigation::update_candidate(navigation.direction(),
                                              navigation.get_target(), track,
                                              det, cfg, ctx)) {
                navigation.set_status(navigation::status::e_unknown);
                navigation.set_fair_trust();
            } else {
                // Update navigation flow on the new candidate information
                navigation::update_status(navigation, cfg);

                navigation.run_inspector(cfg, track.pos(), track.dir(),
                                         "Update complete: high trust: ");

                // The work is done if: the track has not reached a surface yet
                // or trust is gone (portal was reached or the cache is broken).
                if (navigation.get_status() ==
                        navigation::status::e_towards_object ||
                    navigation.get_trust_level() ==
                        navigation::trust_level::e_no_trust) {
                    return false;
                }

                // Else: Track is on module.
                // Ready the next candidate after the current module
                if (navigation::update_candidate(navigation.get_direction(),
                                                 navigation.get_target(), track,
                                                 det, cfg, ctx)) {
                    return false;
                }

                // If next candidate is not reachable, don't 'return', but
                // escalate the trust level.
                // This will run into the fair trust case below.
                navigation.set_fair_trust();
            }
        }

        // Re-evaluate all currently available candidates and sort again
        // - do this when your navigation state is stale, but not invalid
        if (navigation.get_trust_level() == navigation::trust_level::e_fair) {

            for (auto &candidate : navigation) {
                // Disregard this candidate if it is not reachable
                if (!navigation::update_candidate(navigation.get_direction(),
                                                  candidate, track, det, cfg,
                                                  ctx)) {
                    // Forcefully set dist to numeric max for sorting
                    candidate.path = std::numeric_limits<scalar_t>::max();
                }
            }
            detray::sequential_sort(navigation.begin(), navigation.end());
            // Take the nearest (sorted) candidate first
            navigation.set_next(navigation.begin());
            // Ignore unreachable elements (needed to determine exhaustion)
            navigation.set_last(find_invalid(navigation.get_candidates()));
            // Update navigation flow on the new candidate information
            navigation::update_status(navigation, cfg);

            navigation.run_inspector(cfg, track.pos(), track.dir(),
                                     "Update complete: fair trust: ");

            // If a good candidate was found, return (not init)
            if (!navigation.is_exhausted()) {
                return false;
            }
        }

        // Actor flagged cache as broken (other cases of 'no trust' are
        // handeled after volume switch was checked in 'update()')
        if (navigation.trust_level() == navigation::trust_level::e_no_trust) {
            init(track, navigation, cfg, ctx);
            return true;
        }

        return false;
    }

    /// Helper to evict all unreachable/invalid candidates from the cache:
    /// Finds the first unreachable candidate (has been invalidated during
    /// update) in a sorted (!) cache.
    ///
    /// @param candidates the cache of candidates to be cleaned
    DETRAY_HOST_DEVICE inline auto find_invalid(
        const typename caching_navigator::state::candidate_cache_t &candidates)
        const {
        // Depends on previous invalidation of unreachable candidates!
        auto not_reachable = [](const intersection_type &candidate) {
            return candidate.path == std::numeric_limits<scalar_t>::max();
        };

        return detray::find_if(candidates.begin(), candidates.end(),
                               not_reachable);
    }
};

}  // namespace detray
