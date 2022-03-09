/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/core/detector.hpp"
#include "detray/core/intersection.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tools/intersection_kernel.hpp"
#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

namespace navigation {

/** Navigation status flag */
enum status : int {
    e_abort = -3,
    e_exit = -2,
    e_unknown = -1,
    e_towards_object = 0,  // move towards next object
    e_on_target = 1,       // reached object
};

/** Navigation trust level */
enum trust_level : int {
    e_no_trust = 0,    // re-initialize the volume (e.g. local navigation)
    e_fair_trust = 1,  // update the distance & order of the
                       // (preselected) candidates
    e_high_trust = 3,  // update the distance to the next candidate
    e_full_trust = 4   // don't update anything
};

/** A void inpector that does nothing.
 *
 * Inspectors can be plugged in to understand the
 * the current navigation state information.
 *
 */
struct void_inspector {
    template <typename state_t>
    DETRAY_HOST_DEVICE void operator()(const state_t & /*ignored*/,
                                       const char * /*ignored*/) {}
};

}  // namespace navigation

/** The navigator struct is initialized around a detector object, but is itself
 *  agnostic to the object/primitive type. It only requires a link to the next
 *  navigation volume in every candidate that is computed by intersection from
 *  the objects. A module surface should link back to the volume it is conained
 *  in, while a portal surface links to the next volume in the direction of the
 *  track.
 *
 *  This navigator applies a trust level based update of its candidate
 *  (intersection) cache. The trust level, and with it the appropriate update
 *  policy, must be set by an actor. Otherwise, no update will be performed.
 *
 *  It is set up by an init() call and then follows a sequence of
 *  [- step()]
 *  - update()
 *  calls.
 *
 * The heartbeat indicates, that the navigation is still in a valid state.
 *
 * @tparam detector_t the detector to navigate
 * @tparam inspector_t is a validation inspector that can record information
 *         about the navaigation state at different points of the nav. flow.
 */
template <typename detector_t,
          typename inspector_t = navigation::void_inspector>
class navigator {

    public:
    using inspector_type = inspector_t;
    using detector_type = detector_t;
    using volume_type = typename detector_t::volume_type;
    template <typename T>
    using vector_type = typename detector_t::template vector_type<T>;

    /** A navigation state object used to cache the information of the
     *  current navigation stream. These can be read or set in between
     *  navigation calls.
     *
     *  It requires to have a scalar represenation to be used for a stepper
     **/
    class state {
        friend class navigator;

        public:
        /** Default constructor
         **/
        state() = default;

        /** Constructor with memory resource
         **/
        state(vecmem::memory_resource &resource) : _candidates(&resource) {}

        /** Constructor from candidates vector_view
         **/
        DETRAY_HOST_DEVICE state(vector_type<intersection> candidates)
            : _candidates(candidates) {}

        /** Scalar representation of the navigation state,
         * @returns distance to next
         **/
        DETRAY_HOST_DEVICE
        scalar operator()() const { return _next->path; }

        /** @returns current candidates - const */
        DETRAY_HOST_DEVICE
        inline const auto &candidates() const { return _candidates; }

        /** @returns current candidates */
        DETRAY_HOST_DEVICE
        inline auto &candidates() { return _candidates; }

        /** @returns current object that was reached */
        DETRAY_HOST_DEVICE
        inline auto current() const { return _next - 1; }

        /** @returns next object that we want to reach - const */
        DETRAY_HOST_DEVICE
        inline const auto &next() const { return _next; }

        /** @returns next object that we want to reach */
        DETRAY_HOST_DEVICE
        inline auto &next() { return _next; }

        /** Clear the kernel */
        DETRAY_HOST_DEVICE
        void clear() {
            _candidates.clear();
            _next = _candidates.end();
        }

        /** Call the navigation inspector */
        DETRAY_HOST_DEVICE
        inline auto run_inspector(const char *message) {
            return _inspector(*this, message);
        }

        /** @returns the navigation inspector */
        DETRAY_HOST
        inline auto &inspector() { return _inspector; }

        /** @returns current object the navigator is on (might be invalid if
         * between objects) - const
         */
        DETRAY_HOST_DEVICE
        inline auto current_object() const -> dindex { return _object_index; }

        /** Update current object the navigator is on  */
        DETRAY_HOST_DEVICE
        inline void set_object(dindex obj) { _object_index = obj; }

        /** @returns current navigation status - const */
        DETRAY_HOST_DEVICE
        inline auto status() const -> navigation::status { return _status; }

        /** Set new navigation status */
        DETRAY_HOST_DEVICE
        inline void set_status(navigation::status stat) { _status = stat; }

        /** @returns tolerance to determine if we are on object - const */
        DETRAY_HOST_DEVICE
        inline auto tolerance() const -> scalar { return _on_object_tolerance; }

        /** Adjust the on-object tolerance */
        DETRAY_HOST_DEVICE
        inline void set_tolerance(scalar tol) { _on_object_tolerance = tol; }

        /** @returns navigation trust level - const */
        DETRAY_HOST_DEVICE
        inline auto trust_level() const -> navigation::trust_level {
            return _trust_level;
        }

        /** Update navigation trust level */
        DETRAY_HOST_DEVICE
        inline void set_trust_level(navigation::trust_level tr_lvl) {
            _trust_level = tr_lvl;
        }

        /** Update navigation trust level */
        DETRAY_HOST_DEVICE
        inline void set_no_trust() { _trust_level = navigation::e_no_trust; }

        /** Update navigation trust level */
        DETRAY_HOST_DEVICE
        inline void set_full_trust() {
            _trust_level = navigation::e_full_trust;
        }

        /** Update navigation trust level */
        DETRAY_HOST_DEVICE
        inline void set_high_trust() {
            _trust_level = navigation::e_high_trust;
        }

        /** Update navigation trust level */
        DETRAY_HOST_DEVICE
        inline void set_fair_trust() {
            _trust_level = navigation::e_fair_trust;
        }

        /** @returns current volume (index) - const */
        DETRAY_HOST_DEVICE
        inline auto volume() const -> dindex { return _volume_index; }

        /** Set start/new volume */
        DETRAY_HOST_DEVICE
        inline void set_volume(dindex v) { _volume_index = v; }

        /** Helper method to check if a kernel is exhausted - const */
        DETRAY_HOST_DEVICE
        inline auto is_exhausted() const -> bool {
            return (_next == _candidates.end());
        }

        /** Helper method to check if a kernel is exhausted - const */
        template <typename track_t>
        DETRAY_HOST_DEVICE inline auto is_on_object(
            const track_t &track, const intersection &candidate) const -> bool {
            if ((candidate.path < _on_object_tolerance) and
                (candidate.path > track.overstep_tolerance())) {
                return true;
            }
            return false;
        }

        /** Helper method to set common state values */
        DETRAY_HOST_DEVICE
        inline void set_state(const navigation::status status,
                              const dindex obj_idx,
                              const navigation::trust_level trust_level) {
            this->set_object(obj_idx);
            this->set_trust_level(trust_level);
            this->set_status(status);
        }

        /** Navigation state that cannot be recovered from. Leave the other
         *  data for inspection.
         *
         * @return navigation heartbeat (dead)
         */
        DETRAY_HOST_DEVICE
        inline auto abort() -> bool {
            _status = navigation::e_abort;
            // Don't do anything if aborted
            _trust_level = navigation::e_full_trust;
            run_inspector("Aborted: ");
            return false;
        }

        /** Navigation reaches target or leaves detector world. Stop
         *  navigation.
         *
         * @return navigation heartbeat (dead)
         */
        DETRAY_HOST_DEVICE
        inline auto exit() -> bool {
            _status = navigation::e_exit;
            _trust_level = navigation::e_full_trust;
            run_inspector("Exited: ");
            return false;
        }

        /** Check whether navigation was completed
         *
         * @return navigation status
         */
        DETRAY_HOST_DEVICE
        inline auto is_complete() const -> bool {
            // Normal exit for this navigation?
            return _status == navigation::e_exit;
        }

        private:
        /// Our list of candidates (intersections with object)
        vector_type<intersection> _candidates = {};

        /// The next best candidate
        typename vector_type<intersection>::iterator _next = _candidates.end();

        /// Distance to next - will be cast into a scalar with call operator
        scalar _distance_to_next = std::numeric_limits<scalar>::infinity();

        /// The inspector type of this navigation engine
        inspector_type _inspector;

        /// Index of a object (surface/portal) if is reached, otherwise invalid
        dindex _object_index = dindex_invalid;

        /// The navigation status
        navigation::status _status = navigation::e_unknown;

        /// The on object tolerance - permille
        scalar _on_object_tolerance = 1e-3;

        /// The navigation trust level
        navigation::trust_level _trust_level = navigation::e_no_trust;

        /// Volume we are currently navigating in
        dindex _volume_index = 0;
    };

    DETRAY_HOST_DEVICE
    navigator(const detector_t &d) : _detector(&d) {}

    /// Alias call to update for readbility
    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE bool init(state &navigation,
                                 stepper_state_t &stepping) const {
        // Establish navigation status before the first step
        return update(navigation, stepping);
    }

    /** Navigation status() call which established the current navigation
     *  information.
     *
     * @tparam track_t type of the track (including context)
     *
     * @param navigation [in, out] is the navigation state object
     * @param stepping stepper state that contains the track
     *
     * @return a heartbeat to indicate if the navigation is still alive
     **/
    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE inline bool update(state &navigation,
                                          stepper_state_t &stepping) const {
        bool heartbeat = true;

        if (navigation.candidates().capacity() != 0)
            // navigation.run_inspector("Before update: ");

            // If there is no_trust (e.g. at the beginning of the navigation in
            // a volume), the kernel will be initialized. Otherwise the
            // candidates are re-evaluated based on current trust level
            update_kernel(navigation, stepping);

        // Did we hit a portal? Then switch volume
        heartbeat &= check_volume_switch(navigation);

        // If no trust could be restored for the current state, (local)
        // navigation might be exhausted or we switched volumes:
        // (re-)initialize volume
        if (navigation.trust_level() == navigation::e_no_trust) {
            navigation.clear();
            update_kernel(navigation, stepping);
        }

        // Should never be the case after complete update call
        if (navigation.trust_level() != navigation::e_full_trust) {
            heartbeat &= navigation.abort();
        }

        return heartbeat;
    }

    /** Helper method to intersect all objects of a surface/portal store
     *
     * @tparam stepper_state_t state type of the stepper object
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param stepping stepper state that contains the track
     * @param track the track information
     * @param volume the current tracking volume
     */
    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE inline void initialize_kernel(
        state &navigation, stepper_state_t &stepping,
        const volume_type &volume) const {

        const auto &track = stepping();

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
        set_next(track, navigation, stepping);
        // navigation.run_inspector("Init complete: ");
    }

    /** Helper method to update the next candidate intersection based on
     *  trust level. Will initialize the kernel if there is no trust.
     *
     * @tparam stepper_state_t type of the stepper
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param stepping the track information in the stepper state
     */
    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE inline void update_kernel(
        state &navigation, stepper_state_t &stepping) const {

        // Current candidate is up to date, nothing left to do
        if (navigation.trust_level() == navigation::e_full_trust) {
            return;
        }

        const volume_type &volume =
            _detector->volume_by_index(navigation.volume());
        // Navigation state is broken/not initialized
        if (navigation.trust_level() == navigation::e_no_trust) {
            initialize_kernel(navigation, stepping, volume);
            return;
        }

        const auto &track = stepping();
        // Update only the current candidate, or close neighbors if we miss the
        // current candidate - do this only when the navigation state is still
        // largely consistent (high trust)
        if (navigation.trust_level() == navigation::e_high_trust or
            navigation.candidates().size() == 1) {
            while (not navigation.is_exhausted()) {

                intersection &candidate = *navigation.next();
                update_candidate(track, candidate);

                // This is likely the next target
                if (is_reachable(candidate, track)) {

                    if (navigation.is_on_object(track, candidate)) {
                        navigation.set_state(navigation::e_on_target,
                                             candidate.index,
                                             navigation::e_full_trust);
                        // Release stepping constraints
                        stepping.release_step_size();
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
                        navigation.set_state(navigation::e_towards_object,
                                             dindex_invalid,
                                             navigation::e_full_trust);
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
            // navigation.run_inspector("Update (high trust) no candidate
            // found:");
        }
        // Re-evaluate all currently available candidates
        // - do this when your navigation state is stale, but not invalid
        else if (navigation.trust_level() == navigation::e_fair_trust) {
            for (auto &candidate : navigation.candidates()) {
                update_candidate(track, candidate);
                // Disregard this candidate
                if (not is_reachable(candidate, track)) {
                    // Forcefully set dist to numeric max for sorting
                    invalidate_candidate(candidate);
                }
            }
            // Sort again
            set_next(track, navigation, stepping);
            // No surface in this volume is reachable: Kernel is exhausted
            if (not is_reachable(*navigation.next(), track)) {
                // Re-initialize the volume in the next update call
                navigation.set_no_trust();
            }
            // Call the inspector before returning
            // navigation.run_inspector("Update (fair trust): ");
        }
        // Finally, if no candidate was found, trust could no be restored
        if (navigation.is_exhausted()) {
            navigation.set_no_trust();
        }
    }

    /** Helper method to sort within the kernel and find next object
     *
     * @param navigation [in, out] navigation state that contains the kernel
     */
    template <typename track_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE inline void set_next(track_t &track, state &navigation,
                                            stepper_state_t &stepping) const {

        if (not navigation.candidates().empty()) {

            // Sort distance to next & set navigation status
            detail::sequential_sort(navigation.candidates().begin(),
                                    navigation.candidates().end());
            // Remove unreachable elements
            remove_invalid_candidates(navigation.candidates());
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
                navigation.set_state(navigation::e_on_target,
                                     navigation.current()->index,
                                     navigation::e_full_trust);
                // Call the inspector on this portal crossing, then go to next
                navigation.run_inspector("Skipping direct hit: ");
                if (navigation.is_exhausted() or
                    not is_reachable(*navigation.next(), track)) {
                    navigation.set_no_trust();
                    return;
                }
            }
            // Now, we are on our way to the next candidate
            navigation.set_state(navigation::e_towards_object, dindex_invalid,
                                 navigation::e_full_trust);
        } else {
            navigation.set_state(navigation::e_unknown, dindex_invalid,
                                 navigation::e_no_trust);
        }
        // Release stepping constraints
        stepping.release_step_size();
        // Call the inspector on new status
        // navigation.run_inspector("Set next: ");
    }

    /** Helper method to check and perform a volume switch
     *
     * @param navigation is the navigation state
     */
    DETRAY_HOST_DEVICE
    bool check_volume_switch(state &navigation) const {
        // Check if we need to switch volume index and (re-)initialize
        // Do this when we are on a portal that has a volume link different to
        // the volume we currently navigate in
        if (navigation.status() == navigation::e_on_target and
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

    /** Helper method that updates a single candidate
     *
     * @tparam track_t type of the track (including context)
     *
     * @param track the track information
     */
    template <typename track_t>
    DETRAY_HOST_DEVICE inline void update_candidate(
        const track_t &track, intersection &candidate) const {
        const dindex obj_idx = candidate.index;
        candidate =
            intersect(track, _detector->get_surface(obj_idx),
                      _detector->transform_store(), _detector->mask_store());

        candidate.index = obj_idx;
    }

    /** Helper method that invalidates a candidate by setting its path length to
     *  an unreachable value. It then gets sorted to the back of the candidates
     *  cache.
     *
     * @param candidate the candidate to be invalidated
     */
    DETRAY_HOST_DEVICE
    inline void invalidate_candidate(intersection &candidate) const {
        candidate.path = std::numeric_limits<scalar>::max();
    }

    /** Helper to determine if a candidate has been invalidated
     *
     * @param candidate the candidate to be invalidated
     * @returns true if is reachable by track
     */
    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool is_reachable(intersection &candidate,
                                                track_t &track) const {
        return candidate.status == e_inside and
               candidate.path >= track.overstep_tolerance() and
               candidate.path < std::numeric_limits<scalar>::max();
    }

    /** Helper to evict all unreachable/invalid candidates from candidates cache
     *
     * @param candidates the candidates to be cleaned
     */
    template <typename cache_t>
    DETRAY_HOST_DEVICE inline void remove_invalid_candidates(
        cache_t &candidates) const {
        // TODO: Find more performant solution
        auto not_reachable = [](intersection cand) {
            return cand.path == std::numeric_limits<scalar>::max() or
                   cand.status == e_missed;
        };
        auto first_invalid_cand = detail::find_if(
            candidates.begin(), candidates.end(), not_reachable);
        candidates.erase(first_invalid_cand, candidates.end());
    }

    /// @returns pointer to detector
    DETRAY_HOST_DEVICE
    const detector_t &get_detector() const { return *_detector; }

    /** @return the vecmem jagged vector buffer for surface candidates */
    // TODO: det.get_n_max_objects_per_volume() is way too many for
    // candidates size allocation. With the local navigation, the size can be
    // restricted to much smaller value
    DETRAY_HOST
    vecmem::data::jagged_vector_buffer<intersection> create_candidates_buffer(
        const unsigned int n_tracks,
        vecmem::memory_resource &device_resource) const {

        return vecmem::data::jagged_vector_buffer<intersection>(
            std::vector<std::size_t>(n_tracks, 0),
            std::vector<std::size_t>(n_tracks,
                                     _detector->get_n_max_objects_per_volume()),
            device_resource, _detector->resource());
    }

    private:
    /** the containers for all data */
    const detector_t *_detector;
};

/** @return the vecmem jagged vector buffer for surface candidates */
// TODO: det.get_n_max_objects_per_volume() is way too many for
// candidates size allocation. With the local navigation, the size can be
// restricted to much smaller value
template <typename detector_t>
DETRAY_HOST vecmem::data::jagged_vector_buffer<intersection>
create_candidates_buffer(const detector_t &det, const unsigned int n_tracks,
                         vecmem::memory_resource &device_resource) {
    return vecmem::data::jagged_vector_buffer<intersection>(
        std::vector<std::size_t>(n_tracks, 0),
        std::vector<std::size_t>(n_tracks, det.get_n_max_objects_per_volume()),
        device_resource, det.resource());
}

}  // namespace detray