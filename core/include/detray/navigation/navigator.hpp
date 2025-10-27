/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
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
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/navigation/detail/print_state.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_kernel.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/log.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

namespace navigation {

static constexpr std::size_t default_cache_size{8u};

}  // namespace navigation

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
/// This navigator applies a trust level based update of its candidate
/// (intersection) cache, which is kept in the naviagtor's state. The trust
/// level, and with it the appropriate update policy, must be set by an actor,
/// otherwise no update will be performed.
///
/// The navigation state is set up by an init() call and then follows a
/// sequence of
/// - step()       (stepper)
/// - update()     (navigator)
/// - run_actors() (actor chain)
/// - update()     (navigator)
/// calls, which are handled by the propagator class.
///
/// The navigation heartbeat indicates, that the navigation is still running
/// and in a valid state.
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

        friend class navigator;

        // Allow the filling/updating of candidates
        friend struct intersection_initialize<ray_intersector>;
        friend struct intersection_update<ray_intersector>;

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
        /// Insert a new element @param new_cadidate before position @param pos
        DETRAY_HOST_DEVICE
        constexpr void insert(candidate_itr_t pos,
                              const intersection_type &new_cadidate) {

            // Candidate is too far away to be placed in cache
            if (pos == this->candidates().end()) {
                return;
            }

            assert(detail::is_invalid_value(new_cadidate.volume_link) ||
                   new_cadidate.volume_link <
                       this->detector().volumes().size());

            // Insert the first candidate
            if (this->n_candidates() == 0) {
                this->candidates()[0] = new_cadidate;
                this->last_index(this->last_index() + 1);
                assert(this->next_index() <= this->last_index() + 1);
                assert(static_cast<std::size_t>(this->last_index()) <
                       k_cache_capacity);
                return;
            }

            // Position where to insert the new candidate
            auto idx{static_cast<dist_t>(
                detray::ranges::distance(this->candidates().begin(), pos))};
            assert(idx >= 0);

            // Do not add the same surface (intersection) multiple times
            const auto is_clash_at_pos = [this,
                                          &new_cadidate](std::size_t index) {
                return (this->candidates()[index].sf_desc.barcode() ==
                        new_cadidate.sf_desc.barcode()) &&
                       (math::fabs(this->candidates()[index].path() -
                                   new_cadidate.path()) <= 1e-5f);
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
            this->candidates()[static_cast<std::size_t>(idx)] = new_cadidate;
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

    private:
    /// A functor that fills the navigation candidates vector by intersecting
    /// the surfaces in the volume neighborhood
    struct candidate_search {

        /// Test the volume links
        template <typename track_t>
        DETRAY_HOST_DEVICE void operator()(
            const typename detector_type::surface_type &sf_descr,
            const detector_type &det, const context_type &ctx,
            const track_t &track, state &nav_state,
            const intersection::config &intr_cfg) const {

            const auto sf = geometry::surface{det, sf_descr};

            sf.template visit_mask<intersection_initialize<ray_intersector>>(
                nav_state,
                detail::ray<algebra_type>(
                    track.pos(),
                    static_cast<scalar_type>(nav_state.direction()) *
                        track.dir()),
                sf_descr, det.transform_store(), ctx, intr_cfg,
                nav_state.external_tol());
        }
    };

    public:
    /// @brief Helper method to initialize a volume.
    ///
    /// Calls the volumes accelerator structure for local navigation, then tests
    /// the surfaces for intersection and sorts the reachable candidates to find
    /// the clostest one (next candidate).
    ///
    /// @tparam track_t type of track, needs to provide pos() and dir() methods
    ///
    /// @param track access to the track parameters
    /// @param state the current navigation state
    /// @param cfg the navigation configuration
    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx,
        const bool use_path_tolerance_as_overstep_tolerance = true) const {

        DETRAY_VERBOSE_HOST_DEVICE("Called 'init()'");

        const auto &det = navigation.detector();
        const auto volume = tracking_volume{det, navigation.volume()};

        // Do not resurrect a failed/finished navigation state
        assert(navigation.is_alive());
        assert(!track.is_invalid());

        // Clean up state
        navigation.clear_cache();

        // Search for neighboring surfaces and fill candidates into cache
        intersection::config intr_cfg{cfg.intersection};
        intr_cfg.overstep_tolerance = use_path_tolerance_as_overstep_tolerance
                                          ? -cfg.intersection.path_tolerance
                                          : cfg.intersection.overstep_tolerance;
        volume.template visit_neighborhood<candidate_search>(
            track, cfg, ctx, det, ctx, track, navigation, intr_cfg);

        // Determine overall state of the navigation after updating the cache
        update_status(navigation, cfg);

        // If init was not successful, the propagation setup might be broken
        if (navigation.trust_level() != navigation::trust_level::e_full) {
            // Do not exit if backward navigation starts on the outmost portal
            if (navigation.is_on_portal() &&
                navigation.direction() == navigation::direction::e_backward) {
                navigation.trust_level(
                    detail::is_invalid_value(navigation.current().volume_link)
                        ? navigation::trust_level::e_full
                        : navigation::trust_level::e_no_trust);
            } else if (!navigation.is_on_portal()) {
                DETRAY_VERBOSE_HOST_DEVICE("Unable to initialize state...");
            }
        }

        navigation.run_inspector(cfg, track.pos(), track.dir(),
                                 "Init complete: ");
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
    ///
    /// @returns a heartbeat to indicate if the navigation is still alive
    template <typename track_t>
    DETRAY_HOST_DEVICE
#ifdef __CUDA_ARCH__
        __attribute__((always_inline))
#endif
        inline bool
        update(const track_t &track, state &navigation,
               const navigation::config &cfg, const context_type &ctx = {},
               const bool /*is_before_actor*/ = true) const {

        assert(navigation.is_alive());
        assert(!track.is_invalid());

        // Candidates are re-evaluated based on the current trust level.
        // Should result in 'full trust'
        bool is_init = update_kernel(track, navigation, cfg, ctx);

        // Update was completely successful (most likely case)
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            DETRAY_VERBOSE_HOST_DEVICE("Update complete: dist to next %fmm",
                                       navigation());
            return is_init;
        }
        // Otherwise: did we run into a portal?
        else if (navigation.is_on_portal()) {
            // Navigation reached the end of the detector world
            if (detail::is_invalid_value(navigation.current().volume_link)) {
                navigation.exit();
                return is_init;
            }

            // Set volume index to the next volume provided by the portal
            navigation.set_volume(navigation.current().volume_link);

            // Either end of world or valid volume index
            assert(detail::is_invalid_value(navigation.volume()) ||
                   navigation.volume() <
                       navigation.detector().volumes().size());

            // Run inspection when needed (keep for debugging)
            // navigation.run_inspector(cfg, track.pos(), track.dir(), "Volume
            // switch: ");

            init(track, navigation, cfg, ctx);
            is_init = true;

            // Fresh initialization, reset trust and hearbeat even though we are
            // on inner portal
            navigation.trust_level(navigation::trust_level::e_full);

            DETRAY_VERBOSE_HOST_DEVICE("Switched to volume %d",
                                       navigation.volume());
        }
        // If no trust could be restored for the current state, (local)
        // navigation might be exhausted: re-initialize volume
        if (navigation.trust_level() != navigation::trust_level::e_full ||
            navigation.cache_exhausted()) {
            DETRAY_VERBOSE_HOST_DEVICE(
                "Full trust could not be restored! RESCURE MODE: Run init with "
                "large tolerances");

            // Use overstep tolerance instead of path tolerance
            const bool use_path_tolerance_as_overstep_tolerance = false;

            init(track, navigation, cfg, ctx,
                 use_path_tolerance_as_overstep_tolerance);
            is_init = true;

            // Sanity check: Should never be the case after complete update call
            if (navigation.trust_level() != navigation::trust_level::e_full) {
                // Try to save the navigation flow: Look further behind the
                // track
                auto loose_cfg{cfg};
                // Use the max mask tolerance in case a track leaves the volume
                // when a sf is 'sticking' out of the portals due to the tol
                loose_cfg.intersection.overstep_tolerance =
                    math::min(100.f * cfg.intersection.overstep_tolerance,
                              -10.f * cfg.intersection.max_mask_tolerance);

                init(track, navigation, loose_cfg, ctx,
                     use_path_tolerance_as_overstep_tolerance);
            }
        }
        // Unrecoverable
        if (navigation.trust_level() != navigation::trust_level::e_full ||
            navigation.cache_exhausted()) {
            navigation.abort("No reachable surfaces");
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
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update_kernel(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx) const {

        const auto &det = navigation.detector();

        // Current candidates are up to date, nothing left to do
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            return false;
        }

        // Update only the current candidate and the corresponding next target
        // - do this only when the navigation state is still coherent
        if (navigation.trust_level() == navigation::trust_level::e_high) {

            DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - high trust");

            // Update next candidate: If not reachable, 'high trust' is broken
            if (!update_candidate(navigation.direction(), navigation.target(),
                                  track, det, cfg.intersection,
                                  navigation.external_tol(), ctx)) {
                navigation.status(navigation::status::e_unknown);
                navigation.set_fair_trust();
            } else {

                // Update navigation flow on the new candidate information
                update_status(navigation, cfg);

                navigation.run_inspector(cfg, track.pos(), track.dir(),
                                         "Update complete: high trust: ");

                // The work is done if: the track has not reached a surface yet
                // or trust is gone (portal was reached or the cache is broken).
                if (navigation.status() ==
                        navigation::status::e_towards_object ||
                    navigation.trust_level() ==
                        navigation::trust_level::e_no_trust) {
                    return false;
                }

                // Else: Track is on module.
                // Ready the next candidate after the current module
                if (update_candidate(
                        navigation.direction(), navigation.target(), track, det,
                        cfg.intersection, navigation.external_tol(), ctx)) {
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
        if (navigation.trust_level() == navigation::trust_level::e_fair &&
            !navigation.cache_exhausted()) {

            DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - fair trust");

            for (auto &candidate : navigation) {
                // Disregard this candidate if it is not reachable
                if (!update_candidate(navigation.direction(), candidate, track,
                                      det, cfg.intersection,
                                      navigation.external_tol(), ctx)) {
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
            update_status(navigation, cfg);

            navigation.run_inspector(cfg, track.pos(), track.dir(),
                                     "Update complete: fair trust: ");

            if (!navigation.cache_exhausted()) {
                return false;
            }
        }

        // Actor flagged cache as broken (other cases of 'no trust' are
        // handeled after volume switch was checked in 'update()')
        if (navigation.trust_level() == navigation::trust_level::e_no_trust) {

            DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - no trust");

            // Use overstep tolerance instead of path tolerance
            const bool use_path_tolerance_as_overstep_tolerance = false;

            init(track, navigation, cfg, ctx,
                 use_path_tolerance_as_overstep_tolerance);
            return true;
        }

        return false;
    }

    /// @brief Helper method that re-establishes the navigation state after an
    /// update.
    ///
    /// It checks whether the track has reached a surface or is still moving
    /// towards the next surface candidate. If no new next candidate can be
    //  found, it flags 'no trust' in order to trigger a volume initialization.
    ///
    /// @param state the current navigation state
    /// @param cfg the navigation configuration
    DETRAY_HOST_DEVICE inline void update_status(
        state &navigation, const navigation::config &cfg) const {

        // Check whether the track reached the current candidate. Might be a
        // portal, in which case the navigation needs to be re-initialized
        if (!navigation.cache_exhausted() &&
            navigation.has_reached_candidate(navigation.target(), cfg)) {
            navigation.status((navigation.target().sf_desc.is_portal())
                                  ? navigation::status::e_on_portal
                                  : navigation::status::e_on_object);
            // Set the next object that we want to reach (this function is only
            // called once the cache has been updated to a full trust state).
            // Might lead to exhausted cache.
            navigation.advance();
        } else {
            // Otherwise we are moving towards the target
            navigation.status(navigation::status::e_towards_object);
        }
        // Exhaustion happens when after an update no next candidate in the
        // cache is reachable anymore -> triggers init of [new] volume
        // In backwards navigation or with strongly bent tracks, the cache may
        // not be exhausted when trying to exit the volume (the ray is seeing
        // the opposite side of the volume)
        navigation.trust_level(navigation.cache_exhausted() ||
                                       navigation.is_on_portal()
                                   ? navigation::trust_level::e_no_trust
                                   : navigation::trust_level::e_full);
    }

    /// @brief Helper method that updates the intersection of a single candidate
    /// and checks reachability
    ///
    /// @tparam track_t type of track, needs to provide pos() and dir() methods
    ///
    /// @param candidate the candidate intersection to be updated
    /// @param track access to the track parameters
    /// @param det access to the detector (geometry)
    /// @param cfg the navigation configuration
    ///
    /// @returns whether the track can reach this candidate.
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update_candidate(
        const navigation::direction nav_dir, intersection_type &candidate,
        const track_t &track, const detector_type &det,
        const intersection::config &intr_cfg,
        const scalar_type external_mask_tolerance,
        const context_type &ctx) const {

        if (candidate.sf_desc.barcode().is_invalid()) {
            return false;
        }

        const auto sf = geometry::surface{det, candidate.sf_desc};

        // Check whether this candidate is reachable by the track
        return sf.template visit_mask<intersection_update<ray_intersector>>(
            detail::ray<algebra_type>(
                track.pos(), static_cast<scalar_type>(nav_dir) * track.dir()),
            candidate, det.transform_store(), ctx, intr_cfg,
            external_mask_tolerance);
    }

    /// Helper to evict all unreachable/invalid candidates from the cache:
    /// Finds the first unreachable candidate (has been invalidated during
    /// update) in a sorted (!) cache.
    ///
    /// @param candidates the cache of candidates to be cleaned
    DETRAY_HOST_DEVICE inline auto find_invalid(
        typename state::candidate_cache_t &candidates) const {
        // Depends on previous invalidation of unreachable candidates!
        auto not_reachable = [](const intersection_type &candidate) {
            return candidate.path() == std::numeric_limits<scalar_type>::max();
        };

        return detray::find_if(candidates.begin(), candidates.end(),
                               not_reachable);
    }
};

}  // namespace detray
