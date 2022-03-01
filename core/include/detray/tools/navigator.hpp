/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <string>
#include <type_traits>

#include "detray/core/detector.hpp"
#include "detray/core/intersection.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/tools/intersection_kernel.hpp"
#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

/** A void inpector that does nothing.
 *
 * Inspectors can be plugged in to understand the
 * the current navigation state information.
 *
 */
struct void_navigator_inspector {
    template <typename state_t>
    DETRAY_HOST_DEVICE void operator()(const state_t & /*ignored*/,
                                       const char * /*ignored*/) {}
};

/** The navigator struct is agnostic to the object/primitive type. It
 *  only requires a link to the next navigation volume in every candidate
 *  that is computed by intersection from the objects. A module surface should
 *  link back to the volume it is conained in, while a portal surface links
 *  to the next volume in the direction of the track
 *
 * It follows the structure of the Acts::Navigator:
 * a sequence of
 * - status()
 * - target()
 * [- step()]
 * calls.
 *
 * The heartbeat indicates, that the navigation is still in a valid state.
 *
 * @tparam detector_t the detector to navigate
 * @tparam inspector_t is a validation inspector
 */
template <typename detector_t, typename inspector_t = void_navigator_inspector>
class navigator {

    public:
    using detector_type = detector_t;
    using volume_type = typename detector_t::volume_type;
    using inspector_type = inspector_t;

    template <typename T>
    using vector_type = typename detector_t::template vector_type<T>;

    /** Navigation status flag */
    enum navigation_status : int {
        e_on_target = -3,
        e_abort = -2,
        e_unknown = -1,
        e_towards_object = 0,  // move towards next object
        e_on_object = 1,       // reached object
    };

    /** Navigation trust level */
    enum navigation_trust_level : int {
        e_no_trust = 0,    // re-evalute the candidates all over
        e_fair_trust = 1,  // re-evaluate the distance & order of the
                           // (preselected) candidates
        e_high_trust = 3,  // re-evaluate the distance to the next candidate
        e_full_trust = 4   // trust fully: Don't re-evaluate
    };

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
        scalar operator()() const { return _distance_to_next; }

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

        /** Update the distance to next candidate */
        DETRAY_HOST_DEVICE
        inline void set_dist(scalar dist) { _distance_to_next = dist; }

        /** Call the navigation inspector */
        DETRAY_HOST_DEVICE
        inline auto run_inspector(const char *message) {
            return _inspector(*this, message);
        }

        /** @returns the navigation inspector */
        DETRAY_HOST_DEVICE
        inline auto &inspector() { return _inspector; }

        /** @returns current object the navigator is on (might be invalid if
         * between objects) - const
         */
        DETRAY_HOST_DEVICE
        inline auto on_object() const { return _object_index; }

        /** Update current object the navigator is on  */
        DETRAY_HOST_DEVICE
        inline void set_object(dindex obj) { _object_index = obj; }

        /** @returns current navigation status - const */
        DETRAY_HOST_DEVICE
        inline auto status() const { return _status; }

        /** Set new navigation status */
        DETRAY_HOST_DEVICE
        inline void set_status(navigation_status stat) { _status = stat; }

        /** @returns tolerance to determine if we are on object - const */
        DETRAY_HOST_DEVICE
        inline auto tolerance() const { return _on_object_tolerance; }

        /** Adjust the on-object tolerance */
        DETRAY_HOST_DEVICE
        inline void set_tolerance(scalar tol) { _on_object_tolerance = tol; }

        /** @returns navigation trust level - const */
        DETRAY_HOST_DEVICE
        inline auto trust_level() const { return _trust_level; }

        /** Update navigation trust level */
        DETRAY_HOST_DEVICE
        inline void set_trust_level(navigation_trust_level lvl) {
            _trust_level = lvl;
        }

        /** @returns current volume (index) - const */
        DETRAY_HOST_DEVICE
        inline auto volume() const { return _volume_index; }

        /** Set start/new volume */
        DETRAY_HOST_DEVICE
        inline void set_volume(dindex v) { _volume_index = v; }

        /** Helper method to check if a kernel is exhausted - const */
        DETRAY_HOST_DEVICE
        bool is_exhausted() const { return (_next == _candidates.end()); }

        /** Navigation state that cannot be recovered from. Leave the other
         *  data for inspection.
         *
         * @return navigation heartbeat (dead)
         */
        DETRAY_HOST_DEVICE
        inline bool abort() {
            _status = e_abort;
            _trust_level = e_no_trust;
            run_inspector("Aborted: ");
            return false;
        }

        /** Navigation reaches target or leaves detector world. Stop
         *  navigation.
         *
         * @return navigation heartbeat (dead)
         */
        DETRAY_HOST_DEVICE
        inline bool exit() {
            _status = e_on_target;
            _trust_level = e_full_trust;
            run_inspector("Exited: ");
            return false;
        }

        private:
        // Our list of candidates (intersections with object)
        vector_type<intersection> _candidates = {};

        // The next best candidate
        typename vector_type<intersection>::iterator _next = _candidates.end();

        /**  Distance to next - will be cast into a scalar with call operator
         */
        scalar _distance_to_next = std::numeric_limits<scalar>::infinity();

        /** The inspector type of this navigation engine */
        inspector_type _inspector = {};

        /** Index of a object (surface/portal) if is reached, otherwise invalid
         */
        dindex _object_index = dindex_invalid;

        /**  The navigation status */
        navigation_status _status = e_unknown;

        /** The on object tolerance - permille */
        scalar _on_object_tolerance = 1e-3;

        /** The navigation trust level */
        navigation_trust_level _trust_level = e_no_trust;

        /** Volume we are currently navigating in */
        dindex _volume_index = dindex_invalid;
    };

    DETRAY_HOST_DEVICE
    navigator(const detector_t &d) : _detector(&d) {}

    template <typename track_t>
    DETRAY_HOST_DEVICE inline intersection update_next_candidate(
        state &navigation, const track_t &track) const {

        dindex obj_idx = navigation.next()->index;
        auto sfi =
            intersect(track, _detector->surfaces()[obj_idx],
                      _detector->transform_store(), _detector->mask_store());
        sfi.index = obj_idx;
        (*navigation.next()) = sfi;

        return sfi;
    }

    /** Navigation status() call which established the current navigation
     *  information.
     *
     * @tparam stepper_state_t type of stepper state
     *
     * @param navigation [in, out] is the navigation state object
     * @param stepping the stepper state
     *
     * @return a heartbeat to indicate if the navigation is still alive
     **/
    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE inline bool status(state &navigation,
                                          stepper_state_t &stepping) const {
        bool heartbeat = true;

        const auto &track = stepping();

        if (navigation.trust_level() >= e_high_trust) {
            while (not navigation.is_exhausted()) {

                // Only update the next candidate
                auto sfi = update_next_candidate(navigation, track);
                navigation.set_dist(sfi.path);

                if (sfi.status == e_inside) {
                    // We may be on next object (trust level is high)
                    if (std::abs(sfi.path) < navigation.tolerance()) {

                        navigation.set_object(sfi.index);
                        navigation.set_status(e_on_object);
                        navigation.set_trust_level(e_high_trust);

                        ++navigation.next();
                    }
                    // we are certainly not on the next object. Trust fully
                    else {
                        navigation.set_object(dindex_invalid);
                        navigation.set_status(e_towards_object);
                        navigation.set_trust_level(e_full_trust);
                    }

                    // Call the inspector before returning
                    navigation.run_inspector("Update (high trust): ");

                    // Don't sort again when coming from high trust
                    break;
                }
                // If not inside: increase and switch to next
                ++navigation.next();
            }
        }

        // Should never be the case after update call (without portals we are
        // trapped)
        if (navigation.candidates().empty() and heartbeat) {
            return navigation.abort();
        }

        // Did we hit a portal? (kernel needs to be re-initialized next time)
        heartbeat = check_volume_switch(navigation);

        return heartbeat;
    }

    /** Target function of the navigator, finds the next candidates
     *  and set the distance to next
     *
     * @tparam stepper_state_t type of stepper state
     *
     * @param navigation is the navigation state
     * @param stepping the stepper state
     *
     * @return a heartbeat to indicate if the navigation is still alive
     **/
    template <typename stepper_state_t>
    DETRAY_HOST_DEVICE bool target(state &navigation,
                                   stepper_state_t &stepping) const {

        bool heartbeat = true;
        const auto &track = stepping();

        // get surface candidates when there is no trust
        if (navigation.trust_level() == e_no_trust) {

            initialize_kernel(navigation, track,
                              _detector->volumes()[navigation.volume()]);

            stepping.release_step_size(navigation());

            return heartbeat;
        }

        // Loop over all candidates and intersect again all candidates
        // - do this when your trust level is low
        else if (navigation.trust_level() == e_fair_trust) {

            while (not navigation.is_exhausted()) {
                update_next_candidate(navigation, track);
                ++navigation.next();
            }

            set_next(navigation);

            return heartbeat;
        }

        // set the distance to the target surface
        else if (navigation.trust_level() == e_high_trust) {

            while (not navigation.is_exhausted()) {

                // Only update the next candidate
                auto sfi = update_next_candidate(navigation, track);
                navigation.set_dist(sfi.path);

                if (sfi.status == e_inside) {

                    // Navigator goes towards the next target
                    navigation.set_object(dindex_invalid);
                    navigation.set_status(e_towards_object);
                    navigation.set_trust_level(e_full_trust);

                    stepping.release_step_size(navigation());

                    return heartbeat;
                }

                // If not inside: increase and switch to next
                ++navigation.next();
            }
        }

        else if (navigation.trust_level() == e_full_trust) {
            // Do nothing for full trust
            return heartbeat;
        }

        return heartbeat;
    }

    /** Helper method to intersect all objects of a surface/portal store
     *
     * @tparam stepper_state_t type of stepper state
     * @tparam range_t the type of range in the detector data containers
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param stepping the stepper state
     * @param volume the current tracking volume
     *
     */
    template <typename track_t>
    DETRAY_HOST_DEVICE inline void initialize_kernel(
        state &navigation, track_t &track, const volume_type &volume) const {

        detail::call_reserve(navigation.candidates(), volume.n_objects());

        // Loop over all indexed objects in volume, intersect and fill
        // @todo - will come from the local object finder
        for (const auto [obj_idx, obj] :
             enumerate(_detector->surfaces(), volume)) {
            // Retrieve candidate from the object
            auto sfi = intersect(track, obj, _detector->transform_store(),
                                 _detector->mask_store());

            // Candidate is invalid if it oversteps too far (this is neg!)
            if (sfi.path < track.overstep_tolerance()) {
                continue;
            }
            // Accept if inside
            if (sfi.status == e_inside) {
                // object the candidate belongs to
                sfi.index = obj_idx;
                navigation.candidates().push_back(sfi);
            }
        }
        // What is the next object we want to reach?
        set_next(navigation);
        navigation.run_inspector("Init: ");
    }

    /** Helper method to sort within the kernel and find next object
     *
     * @param navigation [in, out] navigation state that contains the kernel
     */
    DETRAY_HOST_DEVICE
    inline void set_next(state &navigation) const {

        // Sort distance to next & set navigation status
        if (not navigation.candidates().empty()) {

            // Take the nearest candidate first
            detail::sequential_sort(navigation.candidates().begin(),
                                    navigation.candidates().end());

            navigation.next() = navigation.candidates().begin();

            // Are we still on an object from a previous navigation pass? Then
            // goto the next candidate.
            // This also excludes adjacent portals -> we are on the next portal
            if (navigation() < navigation.tolerance()) {
                // Set it briefly so that the inspector can catch this state
                navigation.set_object(navigation.next()->index);
                // The next object that we want to reach
                ++navigation.next();
                navigation.set_status(e_on_object);
                // Call the inspector on this portal crossing, then go to next
                navigation.run_inspector("Skipping direct hit: ");
            }

            navigation.set_dist(navigation.next()->path);
            // Generally, we are on our way to some candidate
            navigation.set_status(e_towards_object);
            navigation.set_object(dindex_invalid);
            // This is only called after full (re-)evaluation
            navigation.set_trust_level(e_full_trust);

            // Call the inspector on new status
            navigation.run_inspector("Set next: ");
        }
    }

    /** Helper method to check and perform a volume switch
     *
     * @param navigation is the navigation state
     */
    DETRAY_HOST_DEVICE
    bool check_volume_switch(state &navigation) const {
        // Check if we need to switch volume index and (re-)initialize

        if (navigation.status() == e_on_object and
            navigation.volume() != navigation.current()->link) {

            // Set volume index to the next volume provided by the object
            navigation.set_volume(navigation.current()->link);

            navigation.clear();
            navigation.set_trust_level(e_no_trust);

            // We reached the end of the detector world
            if (navigation.volume() == dindex_invalid) {
                // heartbeat
                return navigation.exit();
            }
        }
        // heartbeat
        return true;
    }

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

}  // namespace detray