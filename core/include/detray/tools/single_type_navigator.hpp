/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>

#include "detray/core/intersection.hpp"
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
struct void_inspector {
    template <typename state_type>
    void operator()(const state_type & /*ignored*/) {}
};

/** The navigator struct that is agnostic to the object/primitive type. It
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
 * @tparam volume_container provides the volumes
 * @tparam object_container provides portals and module surfaces (objects)
 * @tparam transform_container provides the object transforms
 * @tparam mask_container provides the object masks
 * @tparam inspector_type is a validation inspector
 */
template <typename volume_container, typename object_container,
          typename transform_container, typename mask_container,
          typename inspector_type = void_inspector>
class single_type_navigator {

    public:
    using object_t = typename object_container::value_type;

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

    /** A nested navigation kernel struct that holds the current candiates.
     **/
    template <template <typename...> class vector_type = dvector>
    struct navigation_kernel {

        // Our list of candidates (intersections with object)
        vector_type<intersection> candidates = {};

        // The next best candidate
        typename vector_type<intersection>::iterator next = candidates.end();

        /** Indicate that the kernel is empty */
        bool empty() const { return candidates.empty(); }

        /** Forward the kernel size */
        size_t size() const { return candidates.size(); }

        /** Clear the kernel */
        void clear() {
            candidates.clear();
            next = candidates.end();
        }
    };

    /** A navigation state object used to cache the information of the
     *  current navigation stream. These can be read or set in between
     *  navigation calls.
     *
     *  It requires to have a scalar represenation to be used for a stepper
     **/
    class state {
        friend class single_type_navigator;

        public:
        /** Scalar representation of the navigation state,
         * @returns distance to next
         **/
        scalar operator()() const { return _distance_to_next; }

        /** @returns current candidates */
        inline const auto &candidates() const { return _kernel.candidates; }

        /** @returns current candidates */
        inline auto &candidates() { return _kernel.candidates; }

        /** @returns current object that was reached */
        inline decltype(auto) current() { return _kernel.next - 1; }

        /** @returns next object that we want to reach */
        inline auto &next() { return _kernel.next; }

        /** @returns the navigation kernel that contains the candidates */
        inline const auto &kernel() { return _kernel; }

        /** Clear the current kernel */
        inline void clear() { _kernel.clear(); }

        /** Update the distance to next candidate */
        inline void set_dist(scalar dist) { _distance_to_next = dist; }

        /** Call the navigation inspector */
        inline decltype(auto) inspector() { return _inspector(*this); }

        /** @returns current object the navigator is on (might be invalid if
         * between objects)
         */
        inline const auto &on_object() { return _object_index; }

        /** Update current object the navigator is on  */
        inline void set_object(dindex obj) { _object_index = obj; }

        /** @returns current navigation status */
        inline const auto &status() { return _status; }

        /** Set new navigation status */
        inline void set_status(navigation_status stat) { _status = stat; }

        /** @returns tolerance to determine if we are on object */
        inline const auto &tolerance() { return _on_object_tolerance; }

        /** Adjust the on-object tolerance */
        inline void set_tolerance(scalar tol) { _on_object_tolerance = tol; }

        /** @returns navigation trust level */
        inline const auto &trust_level() { return _trust_level; }

        /** Update navigation trust level */
        inline void set_trust_level(navigation_trust_level lvl) {
            _trust_level = lvl;
        }

        /** @returns current volume (index) */
        inline const auto &volume() { return _volume_index; }

        /** Set start/new volume */
        inline void set_volume(dindex v) { _volume_index = v; }

        /** Navigation state that cannot be recovered from. Leave the other
         *  data for inspection.
         *
         * @return navigation heartbeat (dead)
         */
        inline bool abort() {
            _status = e_abort;
            _trust_level = e_no_trust;
            return false;
        }

        /** Navigation reaches target or leaves detector world. Stop
         *  navigation.
         *
         * @return navigation heartbeat (dead)
         */
        inline bool exit() {
            _status = e_on_target;
            _trust_level = e_full_trust;
            return false;
        }

        private:
        /**  Distance to next - will be cast into a scalar with call operator
         */
        scalar _distance_to_next = std::numeric_limits<scalar>::infinity();

        /** Kernel for the objects */
        navigation_kernel<> _kernel;

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

    /** Constructor from collection of data containers
     *
     * @param volumes the container for naviagtion volumes
     * @param objects the container for surfaces/portals
     * @param transforms the container for surface/portal transforms
     * @param masks the container for urface/portal masks ("unrollable")
     */
    single_type_navigator(const volume_container &volumes,
                          const object_container &objects,
                          const transform_container &transforms,
                          const mask_container &masks)
        : _volumes(volumes),
          _objects(objects),
          _transforms(transforms),
          _masks(masks) {}

    /** Navigation status() call which established the current navigation
     *  information.
     *
     * @tparam track_t type of the track (including context)
     *
     * @param navigation [in, out] is the navigation state object
     * @param track [in] is the track infromation
     *
     * @return a heartbeat to indicate if the navigation is still alive
     **/
    template <typename track_t>
    inline bool status(state &navigation, const track_t &track) const {

        bool heartbeat = true;

        // If there is no_trust (e.g. at the beginning of the navigation in a
        // volume), the kernel will be initialized. Otherwise the candidates
        // are re-evaluated based on current trust level
        update_kernel(navigation, track, _volumes[navigation.volume()]);

        // Should never be the case after update call (without portals we are
        // trapped)
        if (navigation.kernel().empty()) {
            return navigation.abort();
        }

        // Did we hit a portal? (kernel needs to be re-initialized next time)
        check_volume_switch(navigation);

        // Call the inspector before returning
        navigation.inspector();

        return heartbeat;
    }

    /** Target function of the navigator, finds the next candidates
     *  and set the distance to next
     *
     * @tparam track_t type of the track (including context)
     *
     * @param navigation is the navigation state
     * @param track is the current track information
     *
     * @return a heartbeat to indicate if the navigation is still alive
     **/
    template <typename track_t>
    bool target(state &navigation, const track_t &track) const {

        // We are already on the right track, nothing left to do
        if (navigation.trust_level() == e_full_trust) {
            // heartbeat
            return true;
        }

        // Re-establish navigation status after [external] changes to
        // the navigation state
        return status(navigation, track);
    }

    /** Helper method to intersect all objects of a surface/portal store
     *
     * @tparam track_t type of the track (including context)
     * @tparam range_t the type of range in the detector data containers
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param track the track information
     * @param volume the current tracking volume
     *
     */
    template <typename track_t>
    inline void initialize_kernel(
        state &navigation, const track_t &track,
        const typename volume_container::value_type &volume) const {

        // Get the max number of candidates & run them through the kernel
        navigation.candidates().reserve(volume.n_objects());

        // Loop over all indexed objects in volume, intersect and fill
        // @todo - will come from the local object finder
        for (const auto [obj_idx, obj] : enumerate(_objects, volume)) {

            // Retrieve candidate from the object
            auto sfi = intersect(track, obj, _transforms, _masks);

            // Candidate is invalid if it oversteps too far (this is neg!)
            if (sfi.path < track.overstep_tolerance) {
                continue;
            }
            // Accept if inside, but not the object we are already on
            if (sfi.status == e_inside) {
                // object the candidate belongs to
                sfi.index = obj_idx;
                // the next volume if we encounter the candidate
                sfi.link = std::get<0>(obj.edge());
                navigation.candidates().push_back(sfi);
            }
        }
        // What is the next object we want to reach?
        set_next(navigation);
    }

    /** Helper method to the update the next candidate intersection. Will
     *  initialize candidates if there is no trust in the current navigation
     *  state.
     *
     * @tparam track_t type of the track (including context)
     * @tparam range the type of range in the detector data containers
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param track the track information
     * @param volume the current tracking volume
     *
     * @return A boolean condition if kernel is exhausted or not
     */
    template <typename track_t>
    inline void update_kernel(
        state &navigation, const track_t &track,
        const typename volume_container::value_type &volume) const {

        if (navigation.trust_level() == e_no_trust) {
            initialize_kernel(navigation, track, volume);
            return;
        }
        // Update current candidate, or close neighbors
        // - do this only when you trust level is high
        else if (navigation.trust_level() >= e_high_trust) {
            while (not is_exhausted(navigation.kernel())) {
                // Only update the next candidate
                dindex obj_idx = navigation.next()->index;
                auto sfi =
                    intersect(track, _objects[obj_idx], _transforms, _masks);
                sfi.index = obj_idx;
                sfi.link = std::get<0>(_objects[obj_idx].edge());
                (*navigation.next()) = sfi;
                navigation.set_dist(sfi.path);

                if (sfi.status == e_inside) {
                    // We may be on next object (trust level is high)
                    if (std::abs(sfi.path) < navigation.tolerance()) {
                        navigation.set_object(obj_idx);
                        navigation.set_status(e_on_object);
                        navigation.set_trust_level(e_high_trust);
                        // Set the next object we want to reach might be end()
                        ++navigation.next();
                    }
                    // we are certainly not on the next object. Trust fully
                    else {
                        navigation.set_status(e_towards_object);
                        navigation.set_trust_level(e_full_trust);
                    }
                    // Don't sort again when coming from high trust
                    return;
                }
                // If not inside: increase and switch to next
                ++navigation.next();
            }
        }
        // Loop over all candidates and intersect again all candidates
        // - do this when your trust level is low
        else if (navigation.trust_level() == e_fair_trust) {
            for (auto &candidate : navigation.candidates()) {
                dindex obj_idx = candidate.index;
                auto sfi =
                    intersect(track, _objects[obj_idx], _transforms, _masks);
                candidate = sfi;
                candidate.index = obj_idx;
                candidate.link = std::get<0>(_objects[obj_idx].edge());
            }
            set_next(navigation);
            return;
        }
        // Do nothing on full trust
        return;
    }

    /** Helper method to sort within the kernel and find next object
     *
     * @param navigation [in, out] navigation state that contains the kernel
     */
    inline void set_next(state &navigation) const {

        auto &kernel = navigation._kernel;

        // Sort distance to next & set navigation status
        if (not kernel.candidates.empty()) {

            // Generally, we are on our way to some candidate
            navigation.set_status(e_towards_object);
            // This is only called after full (re-)evaluation
            navigation.set_trust_level(e_full_trust);

            // Take the nearest candidate first
            std::sort(kernel.candidates.begin(), kernel.candidates.end());
            kernel.next = kernel.candidates.begin();

            // Are we still on an object from a previous navigation pass? Then
            // goto the next candidate.
            // This also excludes adjacent portals -> we are towards obj again
            if (navigation() < navigation.tolerance()) {
                // The object we are on
                navigation.set_object(kernel.next->index);
                // The next object that we want to reach
                ++navigation.next();
                // We might be wrong in this assumption, re-evaluate distance
                // to next in the next pass
                navigation.set_trust_level(e_high_trust);
            }
            // No current object
            else {
                navigation.set_object(dindex_invalid);
            }

            navigation.set_dist(kernel.next->path);
            return;
        }
    }

    /** Helper method to check and perform a volume switch
     *
     * @param navigation is the navigation state
     */
    void check_volume_switch(state &navigation) const {
        // Check if we need to switch volume index and (re-)initialize
        if (navigation.status() == e_on_object and
            navigation.volume() != navigation.current()->link) {

            // Set volume index to the next volume provided by the object
            navigation.set_volume(navigation.current()->link);
            navigation.clear();
            navigation.set_trust_level(e_no_trust);

            // We reached the end of the detector world
            if (navigation.volume() == dindex_invalid) {
                navigation.exit();
            }
        }
    }

    /** Helper method to check if a kernel is exhaused
     *
     * @param kernel the kernel to be checked
     *
     * @return true if the kernel is exhaused
     */
    bool is_exhausted(const navigation_kernel<> &kernel) const {
        return (kernel.next == kernel.candidates.end());
    }

    private:
    /** the containers for all data */
    const volume_container &_volumes;
    const object_container &_objects;
    const transform_container &_transforms;
    const mask_container &_masks;
};

}  // namespace detray