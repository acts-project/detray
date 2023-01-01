/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/algorithms.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/surface_finders/neighborhood_kernel.hpp"
#include "detray/utils/ranges.hpp"

// vecmem include(s)
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

namespace navigation {

/// @enum NavigationDirection
/// The navigation direction is always with
/// respect to a given momentum or direction
enum class direction : int { e_backward = -1, e_forward = 1 };

/// Navigation status flags
enum class status {
    e_abort = -3,          ///< error ocurred, propagation will be aborted
    e_on_target = -2,      ///< navigation exited successfully
    e_unknown = -1,        ///< unknown state/not initialized
    e_towards_object = 0,  ///< move towards next object
    e_on_module = 1,       ///< reached module surface
    e_on_portal = 2,       ///< reached portal surface
};

/// Navigation trust levels determine how the candidates chache is updated
enum class trust_level {
    e_no_trust = 0,  ///< re-initialize the volume (i.e. run local navigation)
    e_fair = 1,      ///< update the distance & order of the candidates
    e_high = 3,  ///< update the distance to the next candidate (current target)
    e_full = 4   ///< don't update anything
};

/// A void inpector that does nothing.
///
/// Inspectors can be plugged in to understand the current navigation state.
struct void_inspector {
    template <typename state_t>
    DETRAY_HOST_DEVICE void operator()(const state_t & /*ignored*/,
                                       const char * /*ignored*/) {}
};

}  // namespace navigation

/// The geometry navigation class.
///
/// The navigator is initialized around a detector object, but is itself
/// agnostic to the detectors's object/primitive types.
/// Within a detector volume, the navigatior will perform a local navigation
/// based on the geometry accelerator structure that is provided by the volume.
/// Once the local navigation is resolved, it moves to the next volume by a
/// portal.
/// To this end, it requires a link to the [next] navigation volume in every
/// candidate that is computed by intersection from the detector objects:
/// A module surface must link back to its mothervolume, while a portal surface
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
/// calls, which are handeled by the propagator class.
///
/// The navigation heartbeat indicates, that the navigation is still running
/// and in a valid state.
///
/// @tparam detector_t the detector to navigate
/// @tparam inspector_t is a validation inspector that can record information
///         about the navaigation state at different points of the nav. flow.
template <typename detector_t,
          typename inspector_t = navigation::void_inspector>
class navigator {

    public:
    using intersection_type = line_plane_intersection;
    using inspector_type = inspector_t;
    using detector_type = detector_t;
    using volume_type = typename detector_t::volume_type;
    template <typename T>
    using vector_type = typename detector_t::template vector_type<T>;

    /// A navigation state object used to cache the information of the
    /// current navigation stream.
    ///
    /// The state is passed between navigation calls and is accessible to the
    /// actors in the propagation, for which it defines the public interface
    /// towards the navigation. The navigator is responsible for updating the
    /// elements  in the state's cache with every navigation call, establishing
    /// 'full trust' again.
    class state {
        friend class navigator;

        using candidate_itr_t =
            typename vector_type<intersection_type>::iterator;
        using const_candidate_itr_t =
            typename vector_type<intersection_type>::const_iterator;

        public:
        using detector_type = navigator::detector_type;
        using intersection_t = typename navigator::intersection_type;

        /// Default constructor
        state() = default;

        state(const detector_type &det) : _detector(&det) {}

        /// Constructor with memory resource
        DETRAY_HOST
        state(const detector_type &det, vecmem::memory_resource &resource)
            : _detector(&det), _candidates(&resource) {}

        /// Constructor from candidates vector_view
        DETRAY_HOST_DEVICE state(const detector_type &det,
                                 vector_type<intersection_type> candidates)
            : _detector(&det), _candidates(candidates) {}

        /// @returns a pointer of detector
        DETRAY_HOST_DEVICE
        auto detector() const { return _detector; }

        /// Scalar representation of the navigation state,
        /// @returns distance to next
        DETRAY_HOST_DEVICE
        scalar operator()() const { return _next->path; }

        DETRAY_HOST_DEVICE
        inline void set_unknown() { _status = navigation::status::e_unknown; }

        DETRAY_HOST_DEVICE
        inline void set_on_module() {
            _status = navigation::status::e_on_module;
        }

        /// Updates the iterator position of the last valid candidate
        DETRAY_HOST_DEVICE
        inline void set_next(candidate_itr_t &&new_next) {
            _next = std::move(new_next);
        }

        /// @returns next object that we want to reach (current target)
        DETRAY_HOST_DEVICE
        inline auto next() -> candidate_itr_t & { return _next; }

        /// @returns next object that we want to reach (current target) - const
        DETRAY_HOST_DEVICE
        inline auto next() const -> const const_candidate_itr_t & {
            return _next;
        }

        /// @returns currently cached candidates - const
        DETRAY_HOST_DEVICE
        inline auto candidates() const
            -> const vector_type<intersection_type> & {
            return _candidates;
        }

        /// @returns currently cached candidates
        DETRAY_HOST_DEVICE
        inline auto candidates() -> vector_type<intersection_type> & {
            return _candidates;
        }

        /// @returns numer of currently cached candidates - const
        DETRAY_HOST_DEVICE
        inline auto n_candidates() const ->
            typename std::iterator_traits<candidate_itr_t>::difference_type {
            return std::distance(_candidates.begin(), _last);
        }

        /// @returns current/previous object that was reached
        DETRAY_HOST_DEVICE
        inline auto current() const -> const_candidate_itr_t {
            return _next - 1;
        }

        /*
        /// @returns next object that we want to reach (current target) - const
        DETRAY_HOST_DEVICE
        inline auto next() const -> const const_candidate_itr_t & {
            return _next;
        }
        */

        /// @returns last valid candidate (by position in the cache) - const
        DETRAY_HOST_DEVICE
        inline auto last() const -> const const_candidate_itr_t & {
            return _last;
        }

        /// @returns the navigation inspector
        DETRAY_HOST
        inline auto &inspector() { return _inspector; }

        /// @returns current volume (index) - const
        DETRAY_HOST_DEVICE
        inline auto volume() const -> dindex { return _volume_index; }

        /// Set start/new volume
        DETRAY_HOST_DEVICE
        inline void set_volume(dindex v) { _volume_index = v; }

        /// @returns current object the navigator is on (might be invalid if
        /// between objects) - const
        DETRAY_HOST_DEVICE
        inline auto current_object() const -> dindex { return _object_index; }

        /// @returns the next object the navigator indends to reach
        DETRAY_HOST_DEVICE
        inline auto next_object() const -> dindex { return _next->barcode; }

        /// @returns current navigation status - const
        DETRAY_HOST_DEVICE
        inline auto status() const -> navigation::status { return _status; }

        /// @returns current navigation direction - const
        DETRAY_HOST_DEVICE
        inline auto direction() const -> navigation::direction {
            return _direction;
        }

        /// Set direction
        DETRAY_HOST_DEVICE
        inline void set_direction(const navigation::direction dir) {
            _direction = dir;
        }

        /// @returns tolerance to determine if we are on object - const
        DETRAY_HOST_DEVICE
        inline auto tolerance() const -> scalar { return _on_object_tolerance; }

        /// Adjust the on-object tolerance
        DETRAY_HOST_DEVICE
        inline void set_tolerance(scalar tol) { _on_object_tolerance = tol; }

        /// @returns navigation trust level - const
        DETRAY_HOST_DEVICE
        inline auto trust_level() const -> navigation::trust_level {
            return _trust_level;
        }

        /// Update navigation trust level to no trust
        DETRAY_HOST_DEVICE
        inline void set_no_trust() {
            _trust_level = navigation::trust_level::e_no_trust;
        }

        /// Update navigation trust level to full trust
        DETRAY_HOST_DEVICE
        inline void set_full_trust() {
            _trust_level = _trust_level <= navigation::trust_level::e_full
                               ? _trust_level
                               : navigation::trust_level::e_full;
        }

        /// Update navigation trust level to high trust
        DETRAY_HOST_DEVICE
        inline void set_high_trust() {
            _trust_level = _trust_level <= navigation::trust_level::e_high
                               ? _trust_level
                               : navigation::trust_level::e_high;
        }

        /// Update navigation trust level to fair trust
        DETRAY_HOST_DEVICE
        inline void set_fair_trust() {
            _trust_level = _trust_level <= navigation::trust_level::e_fair
                               ? _trust_level
                               : navigation::trust_level::e_fair;
        }

        /// Helper method to check the track has reached a module surface
        DETRAY_HOST_DEVICE
        inline auto is_on_module() const -> bool {
            return _status == navigation::status::e_on_module;
        }

        /// Helper method to check the track has reached a portal surface
        DETRAY_HOST_DEVICE
        inline auto is_on_portal() const -> bool {
            return _status == navigation::status::e_on_portal;
        }

        /// Helper method to check if a kernel is exhausted - const
        DETRAY_HOST_DEVICE
        inline auto is_exhausted() const -> bool {
            return std::distance(static_cast<const_candidate_itr_t>(_next),
                                 _last) <= 0;
        }

        /// @returns flag that indicates whether navigation was successful
        DETRAY_HOST_DEVICE
        inline auto is_complete() const -> bool {
            // Normal exit for this navigation?
            return _status == navigation::status::e_on_target and !_heartbeat;
        }

        /// Navigation state that cannot be recovered from. Leave the other
        /// data for inspection.
        ///
        /// @return navigation heartbeat (dead)
        DETRAY_HOST_DEVICE
        inline auto abort() -> bool {
            _status = navigation::status::e_abort;
            _heartbeat = false;
            // Don't do anything if aborted
            _trust_level = navigation::trust_level::e_full;
            run_inspector("Aborted: ");
            return _heartbeat;
        }

        /// Navigation reaches target or leaves detector world. Stop
        /// navigation.
        ///
        /// @return navigation heartbeat (dead)
        DETRAY_HOST_DEVICE
        inline auto exit() -> bool {
            _status = navigation::status::e_on_target;
            _heartbeat = false;
            _trust_level = navigation::trust_level::e_full;
            run_inspector("Exited: ");
            this->clear();
            return _heartbeat;
        }

        private:
        /// Helper method to check if a candidate lies on a surface - const
        template <typename track_t>
        DETRAY_HOST_DEVICE inline auto is_on_object(
            const intersection_type &candidate, const track_t &track) const
            -> bool {
            if ((candidate.path < _on_object_tolerance) and
                (candidate.path > track.overstep_tolerance())) {
                return true;
            }
            return false;
        }

        /// Helper to determine if a candidate has been invalidated
        ///
        /// @param candidate the candidate to be invalidated
        /// @returns true if is reachable by track
        template <typename track_t>
        DETRAY_HOST_DEVICE inline auto is_reachable(
            const intersection_type &candidate, track_t &track) const -> bool {
            return candidate.status == intersection::status::e_inside and
                   candidate.path < std::numeric_limits<scalar>::max() and
                   candidate.path >= track.overstep_tolerance();
        }

        /*
        /// @returns next object that we want to reach (current target)
        DETRAY_HOST_DEVICE
        inline auto next() -> candidate_itr_t & { return _next; }
        */
        /*
        /// Updates the iterator position of the last valid candidate
        DETRAY_HOST_DEVICE
        inline void set_next(candidate_itr_t &&new_next) {
            _next = std::move(new_next);
        }
        */

        /// Updates the iterator position of the last valid candidate
        DETRAY_HOST_DEVICE
        inline void set_last(candidate_itr_t &&new_last) {
            _last = std::move(new_last);
        }

        /// Clear the kernel
        DETRAY_HOST_DEVICE
        inline void clear() {
            _candidates.clear();
            _next = _candidates.end();
        }

        /// Call the navigation inspector
        DETRAY_HOST_DEVICE
        inline void run_inspector(const char *message) {
            _inspector(*this, message);
        }

        /// Helper method to set common state values during naviagtion
        ///
        /// @param status current navigation status
        /// @param obj_idx current object (invalid if not on object)
        /// @param trust_lvl current truct level of next target state
        DETRAY_HOST_DEVICE
        inline void set_state(const navigation::status status,
                              const geometry::barcode obj_idx,
                              const navigation::trust_level trust_lvl) {
            _object_index = obj_idx;
            _trust_level = trust_lvl;
            _status = status;
        }

        /// Heartbeat of this navigation flow signals navigation is alive
        bool _heartbeat = false;

        /// Detector pointer
        const detector_type *const _detector;

        /// Our cache of candidates (intersections with any kind of surface)
        vector_type<intersection_type> _candidates = {};

        /// The next best candidate
        candidate_itr_t _next = _candidates.end();

        /// The last reachable candidate
        const_candidate_itr_t _last = _candidates.end();

        /// The inspector type of this navigation engine
        inspector_type _inspector;

        /// Index of an object (module/portal) if is reached, otherwise invalid
        geometry::barcode _object_index = dindex_invalid;

        /// The navigation status
        navigation::status _status = navigation::status::e_unknown;

        /// The navigation direction
        navigation::direction _direction = navigation::direction::e_forward;

        /// The on object tolerance - permille
        scalar _on_object_tolerance = 1e-3;

        /// The navigation trust level determines how this states cache is to
        /// be updated in the current navigation call
        navigation::trust_level _trust_level =
            navigation::trust_level::e_no_trust;

        /// Index in the detector volume container of current navigation volume
        dindex _volume_index = 0;
    };

    /// Helper method to initialize a volume.
    ///
    /// Calls the volumes accelerator structure for local navigation, then tests
    /// the surfaces for intersection and sorts the reachable candidates to find
    /// the clostest one (next candidate).
    ///
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param propagation contains the stepper and navigator states
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline bool init(propagator_state_t &propagation) const {

        state &navigation = propagation._navigation;
        const auto det = navigation.detector();
        const auto &track = propagation._stepping();
        const auto &volume = det->volume_by_index(navigation.volume());

        // Clean up state
        navigation.clear();
        navigation._heartbeat = true;
        // Get the max number of candidates & run them through the kernel
        // detail::call_reserve(navigation.candidates(), volume.n_objects());
        // @TODO: switch to fixed size buffer
        detail::call_reserve(navigation.candidates(), 20u);

        // Search for neighboring surfaces and fill candidates into cache
        fill_candidates(det, volume, track, navigation.candidates());

        // Sort all candidates and pick the closest one
        detail::sequential_sort(navigation.candidates().begin(),
                                navigation.candidates().end());
        navigation.set_next(navigation.candidates().begin());
        // No unreachable candidates in cache after local navigation
        navigation.set_last(navigation.candidates().end());
        // Determine overall state of the navigation after updating the cache
        update_navigation_state(track, propagation);

        // Run inspection when needed
        if constexpr (not std::is_same_v<inspector_t,
                                         navigation::void_inspector>) {
            navigation.run_inspector("Init complete: ");
        }

        return navigation._heartbeat;
    }

    /// Complete update of the nvaigation flow.
    ///
    /// Restores 'full trust' state to the cadidates cache and checks whether
    /// the track stepped onto a portal and a volume switch is due. If so, or
    /// when the previous update according to the given trust level
    /// failed to restore trust, it performs a complete reinitialization of the
    /// navigation.
    ///
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param propagation contains the stepper and navigator states
    ///
    /// @return a heartbeat to indicate if the navigation is still alive
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline bool update(
        propagator_state_t &propagation) const {

        state &navigation = propagation._navigation;

        // Candidates are re-evaluated based on the current trust level
        update_kernel(propagation);

        // Update was completely successful (most likely case)
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            return navigation._heartbeat;
        }
        // Otherwise: did we run into a portal?
        if (navigation.status() == navigation::status::e_on_portal) {
            // Set volume index to the next volume provided by the portal
            navigation.set_volume(navigation.current()->volume_link);

            // Navigation reached the end of the detector world
            if (navigation.volume() == dindex_invalid) {
                navigation.exit();
                return navigation._heartbeat;
            }
            // Run inspection when needed (keep for debugging)
            /*if constexpr (not std::is_same_v<inspector_t,
                                         navigation::void_inspector>) {
                navigation.run_inspector("Volume switch: ");
            }*/
        }
        // If no trust could be restored for the current state, (local)
        // navigation might be exhausted or we switched volumes:
        // re-initialize volume
        navigation._heartbeat &= init(propagation);

        // Sanity check: Should never be the case after complete update call
        if (navigation.trust_level() != navigation::trust_level::e_full or
            navigation.is_exhausted()) {
            navigation.abort();
        }

        return navigation._heartbeat;
    }

    private:
    /// Helper method to update the candidates (surface intersections)
    /// based on an externally provided trust level. Will (re-)initialize the
    /// navigation if there is no trust.
    ///
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param propagation contains the stepper and navigator states
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void update_kernel(
        propagator_state_t &propagation) const {

        state &navigation = propagation._navigation;
        const auto det = navigation.detector();
        const auto &track = propagation._stepping();

        // Current candidates are up to date, nothing left to do
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            return;
        }

        // Update only the current candidate and the corresponding next target
        // - do this only when the navigation state is still coherent
        if (navigation.trust_level() == navigation::trust_level::e_high or
            navigation.n_candidates() == 1) {

            // Update next candidate: If not reachable, 'high trust' is broken
            if (not update_candidate(*navigation.next(), track, det)) {
                navigation.set_state(navigation::status::e_unknown,
                                     dindex_invalid,
                                     navigation::trust_level::e_no_trust);
                return;
            }

            // Update navigation flow on the new candidate information
            update_navigation_state(track, propagation);

            // Run high trust inspection
            if constexpr (not std::is_same_v<inspector_t,
                                             navigation::void_inspector>) {
                navigation.run_inspector("Update complete: high trust: ");
            }

            // The work is done if: the track has not reached a surface yet or
            // trust is gone (portal was reached or the cache is broken).
            if (navigation.status() == navigation::status::e_towards_object or
                navigation.trust_level() ==
                    navigation::trust_level::e_no_trust) {
                return;
            }

            // Else: Track is on module.
            // Ready the next candidate after the current module
            if (update_candidate(*navigation.next(), track, det)) {
                return;
            }

            // If next candidate is not reachable, don't 'return', but
            // escalate the trust level.
            // This will run into the fair trust case below.
            navigation.set_fair_trust();
        }

        // Re-evaluate all currently available candidates
        // - do this when your navigation state is stale, but not invalid
        if (navigation.trust_level() == navigation::trust_level::e_fair) {

            for (auto &candidate : navigation.candidates()) {
                // Disregard this candidate if it is not reachable
                if (not update_candidate(candidate, track, det)) {
                    // Forcefully set dist to numeric max for sorting
                    candidate.path = std::numeric_limits<scalar>::max();
                }
            }
            // Sort again
            detail::sequential_sort(navigation.candidates().begin(),
                                    navigation.candidates().end());
            // Take the nearest candidate first
            navigation.set_next(navigation.candidates().begin());
            // Ignore unreachable elements (needed to determine exhaustion)
            navigation.set_last(find_invalid(navigation.candidates()));
            // Update navigation flow on the new candidate information
            update_navigation_state(track, propagation);

            // Run fair trust inspection
            if constexpr (not std::is_same_v<inspector_t,
                                             navigation::void_inspector>) {
                navigation.run_inspector("Update complete: fair trust: ");
            }
            return;
        }

        // Actor flagged cache as broken (other cases of 'no trust' are
        // handeled after volume switch was checked in 'update()')
        if (navigation.trust_level() == navigation::trust_level::e_no_trust) {
            navigation._heartbeat &= init(propagation);
            return;
        }
    }

    /// Helper method that re-establishes the navigation state after an update.
    ///
    /// It checks wether the track has reached a surface or is still moving
    /// towards the next surface candidate. If no new next candidate can be
    //  found, it flags 'no trust' in order to trigger a volume initialization.
    /// Additionally, all stepper constraints are lifted if a surface is
    /// reached.
    ///
    /// @tparam track_t the type of the track parametrisation
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param track the track that belongs to the current propagation state
    /// @param propagation contains the stepper and navigator states
    template <typename track_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void update_navigation_state(
        track_t &track, propagator_state_t &propagation) const {

        state &navigation = propagation._navigation;

        // Check wether the track reached the current candidate. Might be a
        // portal, in which case the navigation becomes exhausted (the
        // exit-portal is the last reachable surface in every volume)
        if (navigation.is_on_object(*navigation.next(), track)) {
            // Set the next object that we want to reach (this function is only
            // called once the cache has been updated to a full trust state).
            // Might lead to exhausted cache.
            ++navigation.next();
            // Release actor constraints
            propagation._stepping.release_step();
            // Update state accordingly
            navigation.set_state(
                navigation.volume() != navigation.current()->volume_link
                    ? navigation::status::e_on_portal
                    : navigation::status::e_on_module,
                navigation.current()->barcode, navigation::trust_level::e_full);
        } else {
            // Otherwise the track is moving towards a surface
            navigation.set_state(navigation::status::e_towards_object,
                                 dindex_invalid,
                                 navigation::trust_level::e_full);
        }
        // Generally happens when after an update no next candidate in the
        // cache is reachable anymore -> triggers init of [new] volume
        if (navigation.is_exhausted()) {
            navigation.set_no_trust();
        }
    }

    /// Helper method that updates the intersection of a single candidate and
    /// checks reachability
    ///
    /// @tparam track_t type of the track parametrization
    ///
    /// @param candidate the intersection to be updated
    /// @param track the track information
    ///
    /// @returns whether the track can reach this candidate.
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update_candidate(
        intersection_type &candidate, const track_t &track,
        const detector_type *det) const {

        const auto &mask_store = det->mask_store();
        const auto &sf = det->surfaces(candidate.barcode);
        candidate = mask_store.template visit<intersection_update>(
            sf.mask(), detail::ray(track), sf, det->transform_store());

        // Check whether this candidate is reachable by the track
        return candidate.status == intersection::status::e_inside and
               candidate.path >= track.overstep_tolerance();
    }

    /// Helper method that performs the neighborhood lookup
    ///
    /// @tparam track_t type of the track parametrization
    ///
    /// @param track the track information
    ///
    /// @returns an iterable over the detector surface container
    template <int I = static_cast<int>(volume_type::object_id::e_size) - 1,
              typename track_t, typename cache_t>
    DETRAY_HOST_DEVICE inline void fill_candidates(const detector_type *det,
                                                   const volume_type &volume,
                                                   const track_t &track,
                                                   cache_t &candidates) const {
        const auto &surfaces = det->surface_store();
        const auto &link{volume.template link<
            static_cast<typename volume_type::object_id>(I)>()};

        // Only run the query, if object type is contained in volume
        if (detail::get<1>(link) != dindex_invalid) {
            for (const auto &obj : surfaces.template visit<neighborhood_getter>(
                     link, *det, volume, track)) {

                det->mask_store().template visit<intersection_initialize>(
                    obj.mask(), candidates, detail::ray(track), obj,
                    det->transform_store());
            }
        }
        // Check the next surface type
        if constexpr (I > 0) {
            return fill_candidates<I - 1>(det, volume, track, candidates);
        }
    }

    /// Helper to evict all unreachable/invalid candidates from the cache:
    /// Finds the first unreachable candidate (has been invalidated during
    /// update) in a sorted (!) cache.
    ///
    /// @param candidates the cache of candidates to be cleaned
    template <typename cache_t>
    DETRAY_HOST_DEVICE inline auto find_invalid(cache_t &candidates) const {
        // Depends on previous invalidation of unreachable candidates!
        auto not_reachable = [](intersection_type &candidate) {
            return candidate.path == std::numeric_limits<scalar>::max();
        };

        return detail::find_if(candidates.begin(), candidates.end(),
                               not_reachable);
    }
};

/// @return the vecmem jagged vector buffer for surface candidates
// TODO: det.get_n_max_objects_per_volume() is way too many for
// candidates size allocation. With the local navigation, the size can be
// restricted to much smaller value
template <typename detector_t>
DETRAY_HOST vecmem::data::jagged_vector_buffer<line_plane_intersection>
create_candidates_buffer(
    const detector_t &det, const unsigned int n_tracks,
    vecmem::memory_resource &device_resource,
    vecmem::memory_resource *host_access_resource = nullptr) {
    return vecmem::data::jagged_vector_buffer<line_plane_intersection>(
        std::vector<std::size_t>(n_tracks, 0),
        std::vector<std::size_t>(n_tracks, det.get_n_max_objects_per_volume()),
        device_resource, host_access_resource);
}

}  // namespace detray