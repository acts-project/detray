/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/detail/navigation_base_state.hpp"
#include "detray/navigation/detail/navigation_helpers.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/utils/ranges.hpp"

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
/// This navigator always updates the target surface that is closest to the
/// current track position in the current track direction. Once the surface is
/// reached, the current surface will be kept in a second cache slot. If the
/// target cannot be reached anymore, the navigation stream is re-initialized.
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
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t candidate type
template <typename detector_t,
          typename inspector_t = navigation::void_inspector,
          typename intersection_t =
              intersection2D<typename detector_t::surface_type,
                             typename detector_t::algebra_type, false>>
class navigator {

    public:
    using detector_type = detector_t;
    using context_type = detector_type::geometry_context;

    using algebra_type = typename detector_type::algebra_type;
    using scalar_type = dscalar<algebra_type>;
    using point3_type = dpoint3D<algebra_type>;
    using vector3_type = dvector3D<algebra_type>;

    using intersection_type = intersection_t;
    using inspector_type = inspector_t;

    /// @brief Navigation state that contains the current and next candidate
    ///
    /// The cache can hold two candidates: the target, at pos 1, and the
    /// current surface (if any) at pos 0. Once a surface is reached,
    /// the corresponding candidate is moved to the position at index zero
    /// in the cache and the position at index one is filled with the next
    /// target.
    class state : navigation::base_state<detector_type, 2, inspector_t,
                                         intersection_t> {

        friend class navigator;

        // Allow the filling/updating of candidates
        friend struct intersection_initialize<ray_intersector>;
        friend struct intersection_update<ray_intersector>;

        using base_type = navigation::base_state<detector_type, 2, inspector_t,
                                                 intersection_t>;

        public:
        /// Use common methods of contructing a nvaigation state
        using base_type::base_type;

        private:
        /// Set the next surface that we want to reach
        DETRAY_HOST_DEVICE
        constexpr void advance() { /* Index remains the same */
        }

        /// Insert a new element @param new_cadidate if it is closer than the
        /// current next candidate. This is the cache filling scheme that is
        /// called during 'local_navigation'
        DETRAY_HOST_DEVICE constexpr void insert(
            typename base_type::candidate_itr_t,
            const intersection_type &new_cadidate) {

            assert(detail::is_invalid_value(new_cadidate.volume_link) ||
                   new_cadidate.volume_link <
                       this->detector().volumes().size());

            // Insert the first candidate
            if (new_cadidate < this->candidates()[1]) {
                this->candidates()[1] = new_cadidate;
            }
        }

        /// Call the navigation inspector
        DETRAY_HOST_DEVICE inline void run_inspector(
            [[maybe_unused]] const navigation::config &cfg,
            [[maybe_unused]] const point3_type &track_pos,
            [[maybe_unused]] const vector3_type &track_dir,
            [[maybe_unused]] const char *message) {
            if constexpr (!std::is_same_v<inspector_t,
                                          navigation::void_inspector>) {
                m_inspector(*this, cfg, track_pos, track_dir, message);
            }
        }
    };

    public:
    /// @brief Helper method to initialize a volume.
    ///
    /// Calls the volumes accelerator structure for local navigation, then tests
    /// the surfaces for intersection and keeps the clostest one.
    /// The closest candidate is set as 'next candidate' or 'target'.
    ///
    /// @tparam track_t type of track, needs to provide pos() and dir() methods
    ///
    /// @param track access to the track parameters
    /// @param state the current navigation state
    /// @param cfg the navigation configuration
    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(const track_t &track, state &navigation,
                                        const navigation::config &cfg,
                                        const context_type &ctx) const {
        // Run local navigation in the current volume
        navigation::local_navigation(track, navigation, cfg, ctx);

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
    DETRAY_HOST_DEVICE inline bool update(const track_t &track,
                                          state &navigation,
                                          const navigation::config &cfg,
                                          const context_type &ctx = {}) const {
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
    /// @param state the current navigation state
    /// @param cfg the navigation configuration
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update_impl(const track_t &track,
                                               state &navigation,
                                               const navigation::config &cfg,
                                               const context_type &ctx) const {

        const auto &det = navigation.detector();

        // Current candidates are up to date, nothing left to do
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            return false;
        }
    }
};

}  // namespace detray
