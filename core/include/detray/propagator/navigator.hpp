/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/core/detector.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/utils/enumerate.hpp"

namespace detray {

namespace navigation {

/// Navigation status flags
enum class status {
    e_abort = -3,          ///< error ocurred, propagation will be aborted
    e_exit = -2,           ///< navigation exited successfully
    e_unknown = -1,        ///< unknown state/not initialized
    e_towards_object = 0,  ///< move towards next object
    e_on_target = 1,       ///< reached object
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

/// The navigator struct is initialized around a detector object, but is itself
/// agnostic to the object/primitive type. It only requires a link to the next
/// navigation volume in every candidate that is computed by intersection from
/// the objects. A module surface should link back to the volume it is conained
/// in, while a portal surface links to the next volume in the direction of the
/// track.
///
/// This navigator applies a trust level based update of its candidate
/// (intersection) cache. The trust level, and with it the appropriate update
/// policy, must be set by an actor. Otherwise, no update will be performed.
///
/// It is set up by an init() call and then follows a sequence of
/// [- step()]
/// - update()
/// calls.
///
/// The heartbeat indicates, that the navigation is still running and in a
/// valid state.
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
    /// current navigation stream. These can be read or set in between
    /// navigation calls.
    ///
    /// It requires to have a scalar represenation to be used for a stepper
    class state {
        friend class navigator;

        using candidate_itr_t =
            typename vector_type<intersection_type>::iterator;

        public:
        /// Default constructor
        state() = default;

        /// Constructor with memory resource
        DETRAY_HOST
        state(vecmem::memory_resource &resource) : _candidates(&resource) {}

        /// Constructor from candidates vector_view
        DETRAY_HOST_DEVICE state(vector_type<intersection_type> candidates)
            : _candidates(candidates) {}

        /// Scalar representation of the navigation state,
        /// @returns distance to next
        DETRAY_HOST_DEVICE
        scalar operator()() const { return _next->path; }

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

        /// @returns current/previous object that was reached
        DETRAY_HOST_DEVICE
        inline auto current() const -> const candidate_itr_t {
            return _next - 1;
        }

        /// @returns next object that we want to reach (current target) - const
        DETRAY_HOST_DEVICE
        inline auto next() const -> const candidate_itr_t & { return _next; }

        /// @returns next object that we want to reach (current target)
        DETRAY_HOST_DEVICE
        inline auto next() -> candidate_itr_t & { return _next; }

        /// @returns last valid candidate (by position in the cache) - const
        DETRAY_HOST_DEVICE
        inline auto last() const -> const candidate_itr_t & { return _last; }

        /// @returns last valid candidate (by position in the cache)
        DETRAY_HOST_DEVICE
        inline auto last() -> candidate_itr_t & { return _last; }

        /// Updates the iterator position of the last valid candidate
        DETRAY_HOST_DEVICE
        inline void set_last(candidate_itr_t &&new_last) {
            _last = std::move(new_last);
        }

        /// Clear the kernel
        DETRAY_HOST_DEVICE
        void clear() {
            _candidates.clear();
            _next = _candidates.end();
        }

        /// Call the navigation inspector
        DETRAY_HOST_DEVICE
        inline void run_inspector(const char *message) {
            _inspector(*this, message);
        }

        /// @returns the navigation inspector
        DETRAY_HOST
        inline auto &inspector() { return _inspector; }

        /// @returns current object the navigator is on (might be invalid if
        /// between objects) - const
        DETRAY_HOST_DEVICE
        inline auto current_object() const -> dindex { return _object_index; }

        /// Update current object the navigator is on
        DETRAY_HOST_DEVICE
        inline void set_object(dindex obj) { _object_index = obj; }

        /// @returns current navigation status - const
        DETRAY_HOST_DEVICE
        inline auto status() const -> navigation::status { return _status; }

        /// Set new navigation status
        DETRAY_HOST_DEVICE
        inline void set_status(navigation::status stat) { _status = stat; }

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

        /// @returns current volume (index) - const
        DETRAY_HOST_DEVICE
        inline auto volume() const -> dindex { return _volume_index; }

        /// Set start/new volume
        DETRAY_HOST_DEVICE
        inline void set_volume(dindex v) { _volume_index = v; }

        /// Helper method to check if a kernel is exhausted - const
        DETRAY_HOST_DEVICE
        inline auto is_exhausted() const -> bool {
            return (_next == this->last());
        }

        /// Helper method to check if a candidate lies on a surface - const
        template <typename track_t>
        DETRAY_HOST_DEVICE inline auto is_on_object(
            const track_t &track, const intersection_type &candidate) const
            -> bool {
            if ((candidate.path < _on_object_tolerance) and
                (candidate.path > track.overstep_tolerance())) {
                return true;
            }
            return false;
        }

        /// Helper method to set common state values
        ///
        /// @param status current navigation status
        /// @param obj_idx current object (invalid if not on object)
        /// @param trust_lvl current truct level of next target state
        DETRAY_HOST_DEVICE
        inline void set_state(const navigation::status status,
                              const dindex obj_idx,
                              const navigation::trust_level trust_lvl) {
            this->set_object(obj_idx);
            this->set_trust_level(trust_lvl);
            this->set_status(status);
        }

        /// Navigation state that cannot be recovered from. Leave the other
        /// data for inspection.
        ///
        /// @return navigation heartbeat (dead)
        DETRAY_HOST_DEVICE
        inline auto abort() -> bool {
            _status = navigation::status::e_abort;
            // Don't do anything if aborted
            _trust_level = navigation::trust_level::e_full;
            run_inspector("Aborted: ");
            return false;
        }

        /// Navigation reaches target or leaves detector world. Stop
        /// navigation.
        ///
        /// @return navigation heartbeat (dead)
        DETRAY_HOST_DEVICE
        inline auto exit() -> bool {
            _status = navigation::status::e_exit;
            _trust_level = navigation::trust_level::e_full;
            run_inspector("Exited: ");
            this->clear();
            return false;
        }

        /// @returns flag that indicates whether navigation was successful
        DETRAY_HOST_DEVICE
        inline auto is_complete() const -> bool {
            // Normal exit for this navigation?
            return _status == navigation::status::e_exit;
        }

        private:
        /// Update navigation trust level - custom
        DETRAY_HOST_DEVICE
        inline void set_trust_level(navigation::trust_level tr_lvl) {
            _trust_level = tr_lvl;
        }

        /// Our cache of candidates (intersections with any kind of surface)
        vector_type<intersection_type> _candidates = {};

        /// The next best candidate
        candidate_itr_t _next = _candidates.end();

        /// The last reachable candidate
        candidate_itr_t _last = _candidates.end();

        /// The inspector type of this navigation engine
        inspector_type _inspector;

        /// Index of a object (surface/portal) if is reached, otherwise invalid
        dindex _object_index = dindex_invalid;

        /// The navigation status
        navigation::status _status = navigation::status::e_unknown;

        /// The on object tolerance - permille
        scalar _on_object_tolerance = 1e-3;

        /// The navigation trust level determines how this states cache is to
        /// be updated in the current navigation call
        navigation::trust_level _trust_level =
            navigation::trust_level::e_no_trust;

        /// Index in the detector volume container of current navigation volume
        dindex _volume_index = 0;
    };

    /// Constructor from detector object, which is not owned by the navigator
    /// and needs to be guaranteed to have a lifetime beyond that of the
    /// navigator
    DETRAY_HOST_DEVICE
    navigator(const detector_t &d) : _detector(&d) {}

    /// Alias call to update for readbility
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE bool init(propagator_state_t &propagation) const {
        // Establish navigation status before the first step
        return update(propagation);
    }

    /// Updates the navigator state, in particular the intersections in the
    /// internal cache, after an external update. After the cache is updated, it
    /// checks whether the state stepped onto a portal and a volume switch is
    /// due. If so, or the previous update according to the given trust level
    /// failed, it performs a complete reinitialization of the navigation.
    ///
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param propagation contains the stepper and navigator states
    ///
    /// @return a heartbeat to indicate if the navigation is still alive
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline bool update(
        propagator_state_t &propagation) const {
        bool heartbeat = true;
        state &navigation = propagation._navigation;

        // If there is no_trust (e.g. at the beginning of the navigation in
        // a volume), the kernel will be initialized. Otherwise the
        // candidates are re-evaluated based on current trust level
        update_kernel(propagation);

        // Did we hit a portal? Then switch volume
        heartbeat &= check_volume_switch(propagation._navigation);

        // If no trust could be restored for the current state, (local)
        // navigation might be exhausted or we switched volumes:
        // (re-)initialize volume
        if (navigation.trust_level() == navigation::trust_level::e_no_trust) {
            navigation.clear();
            update_kernel(propagation);
        }

        // Should never be the case after complete update call
        if (navigation.trust_level() != navigation::trust_level::e_full) {
            heartbeat &= navigation.abort();
        }

        return heartbeat;
    }

    /// Helper method to initialize a volume by calling the volumes accelerator
    /// structure for local navigation. (Currently brute force)
    ///
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param propagation contains the stepper and navigator states
    /// @param volume the current tracking volume to be initialized
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void initialize_kernel(
        propagator_state_t &propagation, const volume_type &volume) const {

        state &navigation = propagation._navigation;
        const auto &track = propagation._stepping();

        // Get the max number of candidates & run them through the kernel
        detail::call_reserve(navigation.candidates(), volume.n_objects());

        // Loop over all indexed objects in volume, intersect and fill
        // @todo - will come from the local object finder
        for (const auto [obj_idx, obj] :
             enumerate(_detector->surfaces(), volume)) {

            auto candidate = intersect(track, obj, _detector->transform_store(),
                                       _detector->mask_store());
            // Accept potential next target
            if (is_reachable(candidate, track)) {
                candidate.index = obj_idx;
                navigation.candidates().push_back(candidate);
            }
        }
        // What is the closest object we can reach?
        set_next(track, propagation);
        // Call the inspector before returning (only for debugging)
        navigation.run_inspector("Init complete: ");
    }

    /// Helper method to update the candidate cache (surface intersections)
    /// based on an externally provoded trust level. Will (re-)initialize the
    /// navigation if there is no trust.
    ///
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param propagation contains the stepper and navigator states
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void update_kernel(
        propagator_state_t &propagation) const {
        state &navigation = propagation._navigation;

        // Current candidate is up to date, nothing left to do
        if (navigation.trust_level() == navigation::trust_level::e_full) {
            return;
        }

        const volume_type &volume =
            _detector->volume_by_index(navigation.volume());
        // Navigation state is broken/not initialized
        if (navigation.trust_level() == navigation::trust_level::e_no_trust) {
            initialize_kernel(propagation, volume);
            return;
        }

        const auto &track = propagation._stepping();
        // Update only the current candidate, or close neighbors if we miss the
        // current candidate - do this only when the navigation state is still
        // largely consistent (high trust)
        if (navigation.trust_level() == navigation::trust_level::e_high or
            navigation.candidates().size() == 1) {
            while (not navigation.is_exhausted()) {

                intersection_type &candidate = *navigation.next();
                update_candidate(track, candidate);

                // This is likely the next target
                if (is_reachable(candidate, track)) {

                    if (navigation.is_on_object(track, candidate)) {
                        navigation.set_state(navigation::status::e_on_target,
                                             candidate.index,
                                             navigation::trust_level::e_full);
                        // Release stepping constraint from actors
                        propagation._stepping.release_step();
                        // Set the next object we want to reach - might be end()
                        ++navigation.next();
                        if (not navigation.is_exhausted()) {
                            navigation.run_inspector(
                                "Update next (high trust):");
                            // Ready the next candidate in line
                            continue;
                        } else {
                            // There is no next candidate left - re-evaluate
                            navigation.set_no_trust();
                            navigation.run_inspector(
                                "Update (high trust, kernel exhausted):");
                            return;
                        }
                    }
                    // we are not on the next object
                    else {
                        navigation.set_state(
                            navigation::status::e_towards_object,
                            dindex_invalid, navigation::trust_level::e_full);
                        navigation.run_inspector("Update (high trust):");
                    }
                    // Don't sort again when coming from high trust
                    return;
                }
                // Try next candidate
                ++navigation.next();
            }
            // Kernel is exhausted at this point
            // Call the inspector before returning (only for debugging)
            navigation.run_inspector("Update (high trust) no candidate found:");
        }
        // Re-evaluate all currently available candidates
        // - do this when your navigation state is stale, but not invalid
        else if (navigation.trust_level() == navigation::trust_level::e_fair) {
            for (auto &candidate : navigation.candidates()) {
                update_candidate(track, candidate);
                // Disregard this candidate
                if (not is_reachable(candidate, track)) {
                    // Forcefully set dist to numeric max for sorting
                    invalidate_candidate(candidate);
                }
            }
            // Sort again
            set_next(track, propagation);
            // No surface in this volume is reachable: Kernel is exhausted
            if (not is_reachable(*navigation.next(), track)) {
                // Re-initialize the volume in the next update call
                navigation.set_no_trust();
            }
            // Call the inspector before returning (only for debugging)
            navigation.run_inspector("Update (fair trust): ");
        }
        // Finally, if no candidate was found, trust could no be restored
        if (navigation.is_exhausted()) {
            navigation.set_no_trust();
        }
    }

    /// Helper method to sort within the cache and find the object which to
    /// target next (shortest distance by straight line)
    ///
    /// @tparam propagator_state_t state type of the propagator
    ///
    /// @param propagation contains the stepper and navigator states
    template <typename track_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void set_next(
        track_t &track, propagator_state_t &propagation) const {
        state &navigation = propagation._navigation;

        if (not navigation.candidates().empty()) {

            // Sort distance to next & set navigation status
            detail::sequential_sort(navigation.candidates().begin(),
                                    navigation.candidates().end());

            // Re-sort after fair trust can result in many invalidated
            // candidates - ignore them
            if (navigation.trust_level() == navigation::trust_level::e_fair) {
                // Remove unreachable elements
                navigation.set_last(find_invalid(navigation.candidates()));
            } else {
                navigation.set_last(navigation.candidates().end());
            }

            // Take the nearest candidate first
            navigation.next() = navigation.candidates().begin();

            // Are we still on an object from a previous navigation pass? Then
            // goto the next candidate.
            // This mainly updates adjacent portals -> we are automatically on
            // the next portal in the new volume
            if (navigation.is_on_object(track, *navigation.next())) {
                // Set the next object that we want to reach (cache is sorted)
                ++navigation.next();
                // Set temporarily, so that the inspector can catch this state
                navigation.set_state(navigation::status::e_on_target,
                                     navigation.current()->index,
                                     navigation::trust_level::e_full);
                // Release actor constraints
                propagation._stepping.release_step();
                // Call the inspector on this portal crossing, then go to next
                navigation.run_inspector("Skipping direct hit: ");
                if (navigation.is_exhausted() or
                    not is_reachable(*navigation.next(), track)) {
                    navigation.set_no_trust();
                    return;
                }
            }
            // Now, we are on our way to the next candidate
            navigation.set_state(navigation::status::e_towards_object,
                                 dindex_invalid,
                                 navigation::trust_level::e_full);
        } else {
            navigation.set_state(navigation::status::e_unknown, dindex_invalid,
                                 navigation::trust_level::e_no_trust);
        }
        // Release stepping constraints from actors
        // Call the inspector before returning (only for debugging)
        navigation.run_inspector("Set next: ");
    }

    /// Helper method to check and perform a volume switch
    ///
    /// @param navigation is the navigation state
    DETRAY_HOST_DEVICE
    bool check_volume_switch(state &navigation) const {
        // Check if we need to switch volume index and (re-)initialize
        // Do this when we are on a portal that has a volume link different to
        // the volume we currently navigate in
        if (navigation.status() == navigation::status::e_on_target and
            navigation.volume() != navigation.current()->link) {
            // Set volume index to the next volume provided by the portal
            navigation.set_volume(navigation.current()->link);
            // Initialize the new volume in the next update call
            navigation.set_no_trust();

            // We reached the end of the detector world
            if (navigation.volume() == dindex_invalid) {
                // heartbeat
                return navigation.exit();
            }
        }
        return true;
    }

    /// Helper method that updates the intersection of a single candidate
    ///
    /// @tparam track_t type of the track (including context)
    ///
    /// @param track the track information
    template <typename track_t>
    DETRAY_HOST_DEVICE inline void update_candidate(
        const track_t &track, intersection_type &candidate) const {
        const dindex obj_idx = candidate.index;
        candidate =
            intersect(track, _detector->surface_by_index(obj_idx),
                      _detector->transform_store(), _detector->mask_store());

        candidate.index = obj_idx;
    }

    /// Helper method that invalidates a candidate by setting its path length to
    /// an unreachable value. It then gets sorted to the back of the candidates
    /// cache.
    ///
    /// @param candidate the candidate to be invalidated
    DETRAY_HOST_DEVICE
    inline void invalidate_candidate(intersection_type &candidate) const {
        candidate.path = std::numeric_limits<scalar>::max();
    }

    /// Helper to determine if a candidate has been invalidated
    ///
    /// @param candidate the candidate to be invalidated
    /// @returns true if is reachable by track
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool is_reachable(intersection_type &candidate,
                                                track_t &track) const {
        return candidate.status == intersection::status::e_inside and
               candidate.path >= track.overstep_tolerance() and
               candidate.path < std::numeric_limits<scalar>::max();
    }

    /// Helper to evict all unreachable/invalid candidates from the cache
    ///
    /// @param candidates the cache of candidates to be cleaned
    template <typename cache_t>
    DETRAY_HOST_DEVICE inline auto find_invalid(cache_t &candidates) const {
        // Member functions cannot be used here easily (?)
        auto not_reachable = [](intersection_type &candidate) {
            return candidate.path == std::numeric_limits<scalar>::max();
        };

        return detail::find_if(candidates.begin(), candidates.end(),
                               not_reachable);
    }

    /// @returns access to detector
    DETRAY_HOST_DEVICE
    const detector_t &get_detector() const { return *_detector; }

    private:
    /// the containers for all data
    const detector_t *const _detector;
};

/// @return the vecmem jagged vector buffer for surface candidates
// TODO: det.get_n_max_objects_per_volume() is way too many for
// candidates size allocation. With the local navigation, the size can be
// restricted to much smaller value
template <typename detector_t>
DETRAY_HOST vecmem::data::jagged_vector_buffer<line_plane_intersection>
create_candidates_buffer(const detector_t &det, const unsigned int n_tracks,
                         vecmem::memory_resource &device_resource) {
    return vecmem::data::jagged_vector_buffer<line_plane_intersection>(
        std::vector<std::size_t>(n_tracks, 0),
        std::vector<std::size_t>(n_tracks, det.get_n_max_objects_per_volume()),
        device_resource, det.resource());
}

}  // namespace detray