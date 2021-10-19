/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>

#include "core/intersection.hpp"
#include "tools/intersection_kernel.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

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
 *  that is computed by intersection from the objects.
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
 * @tparam data_container is the type of the conatiner that provides the
 *                        geometry data.
 * @tparam inspector_type is a validation inspector
 */
template <typename volume_container, typename object_container,
          typename transform_container, typename mask_container,
          typename inspector_type = void_inspector>
class single_type_navigator {

    public:
    using object_t = typename object_container::value_type;
    using link_t = typename object_t::edge_links;

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
    template <template <typename> class vector_type = dvector>
    struct navigation_kernel {
        // Where are we (nullptr if we are in between objects)
        const object_t *on = nullptr;

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
            on = nullptr;
        }
    };

    /** A navigation state object used to cache the information of the
     * current navigation stream. These can be read or set in between
     * navigation calls.
     *
     * It requires to have a scalar represenation to be used for a stepper
     **/
    class state {
        friend class single_type_navigator;

        public:
        /** Scalar representation of the navigation state,
         * @returns distance to next
         **/
        scalar operator()() const { return distance_to_next; }

        /** Current kernel */
        const auto &nav_kernel() { return kernel; }

        /** Current candidates */
        const auto &candidates() { return kernel.candidates; }

        /** Current volume */
        const auto &volume() { return volume_index; }

        /** Start volume */
        void set_initial_volume(dindex initial_volume) {
            volume_index = initial_volume;
        }

        /** Tolerance to determine of we are on object */
        const auto &tolerance() { return on_object_tolerance; }

        /** Adjust the on-object tolerance */
        void set_tolerance(scalar tol) { on_object_tolerance = tol; }

        /** get the navigation inspector */
        const auto &nav_inspector() { return inspector; }

        /** Current navigation status */
        const auto &nav_status() { return status; }

        /** Current object the navigator is on (might be invalid if between
         * objects)
         */
        const auto &on_object() { return object_index; }

        /** The links (next volume, next object finder) of current
         * candidate
         */
        auto &nav_links() { return links; }

        /** Navigation trust level */
        const auto &nav_trust_level() { return trust_level; }

        /** Navigation trust level */
        void set_trust_level(navigation_trust_level lvl) { trust_level = lvl; }

        const auto &next() { return kernel.next; }
        inline const decltype(auto) current() { return kernel.next - 1; }

        private:
        /** Navigation state cannot be recovered from. Leave the other
         *  data for inspection.
         *
         * @return navigation heartbeat (dead)
         */
        bool abort() {
            status = e_abort;
            trust_level = e_no_trust;
            return false;
        }

        /** Navigation reaches target or leaves detector world. Stop
         *  navigation.
         *
         * @return navigation heartbeat (dead)
         */
        bool exit() {
            status = e_on_target;
            trust_level = e_full_trust;
            return false;
        }

        /** Kernel for the objects */
        navigation_kernel<> kernel;

        /** Volume we are currently navigating in */
        dindex volume_index = dindex_invalid;

        /**  Distance to next - will be cast into a scalar with call operator
         */
        scalar distance_to_next = std::numeric_limits<scalar>::infinity();

        /** The on object tolerance - permille */
        scalar on_object_tolerance = 1e-3;

        /** The inspector type of this navigation engine */
        inspector_type inspector = {};

        /**  The navigation status */
        navigation_status status = e_unknown;

        /** Index of a object (surface/portal) if is reached, otherwise
         * invalid
         */
        dindex object_index = dindex_invalid;

        // Point to the next volume and object finder (not needed here)
        link_t links = {};

        /** The navigation trust level */
        navigation_trust_level trust_level = e_no_trust;
    };

    /** Constructor with move constructor
     *
     * @param data the container for all data: volumes, primitives (objects
     *  i.e. surfaces and portals), transforms and masks (in a tuple by
     *  type)
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
        update_kernel(navigation, track,
                      _volumes[navigation.volume_index].full_range());

        // Should never be the case after update call (without portals we are
        // trapped)
        if (navigation.kernel.empty()) {
            return navigation.abort();
        }
        // Did we hit a portal? (kernel needs to be re-initialized next time)
        check_volume_switch(navigation);
        // Call the inspector before returning
        navigation.inspector(navigation);

        return heartbeat;
    }

    /** Target function of the navigator, finds the next candidates
     *  and set the distance to next
     *
     * @param navigation is the navigation state
     * @param track is the current track information
     *
     * @return a heartbeat to indicate if the navigation is still alive
     **/
    template <typename track_t>
    bool target(state &navigation, const track_t &track) const {

        // We are already on the right track, nothing left to do
        if (navigation.trust_level == e_full_trust) {
            // heartbeat
            return true;
        }

        // Re-establish navigation status after [external] changes to
        // the navigation state
        return status(navigation, track);
    }

    /** Helper method to intersect all objects of a surface/portal store
     *
     * @tparam range the type of range in the detector data containers
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param track the track information
     * @param obj_range the surface/portal index range in the detector cont
     *
     */
    template <typename track_t, typename range_t>
    inline void initialize_kernel(state &navigation, const track_t &track,
                                  const range_t &obj_range) const {

        // Get the max number of candidates & run them through the kernel
        navigation.kernel.candidates.reserve(obj_range[1] - obj_range[0]);

        // Loop over all indexed objects, intersect and fill
        // @todo - will come from the local object finder
        for (size_t obj_idx = obj_range[0]; obj_idx < obj_range[1]; obj_idx++) {
            // Get next object
            const auto &obj = _objects[obj_idx];

            // Retrieve candidate from the object
            auto sfi = intersect(track, obj, _transforms, _masks,
                                 navigation.nav_links());

            // Candidate is invalid if it oversteps too far (this is neg!)
            if (sfi.path < track.overstep_tolerance) {
                continue;
            }
            // Accept if inside, but not the object we are already on
            if (sfi.status == e_inside) {
                // object the candidate belongs to
                sfi.index = obj_idx;
                // the next volume if we encounter the candidate
                sfi.link = obj.edge()[0];
                navigation.kernel.candidates.push_back(sfi);
            }
        }
        // What is the next object we want to reach?
        set_next(navigation);
    }

    /** Helper method to the update the next candidate intersection
     *
     * @tparam range the type of range in the detector data containers
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param track the track information
     * @param obj_range the surface/portal index range in the detector cont
     *
     * @return A boolean condition if kernel is exhausted or not
     */
    template <typename track_t, typename range_t>
    void update_kernel(state &navigation, const track_t &track,
                       const range_t &obj_range) const {

        if (navigation.trust_level == e_no_trust) {
            // This kernel cannot be trusted
            initialize_kernel(navigation, track, obj_range);
            return;
        }
        // Update current candidate, or close neighbors
        // - do this only when you trust level is high
        else if (navigation.trust_level >= e_high_trust) {
            while (not is_exhausted(navigation.kernel)) {
                // Only update the next candidate
                dindex obj_idx = navigation.kernel.next->index;
                const auto &obj = _objects[obj_idx];
                auto sfi = intersect(track, obj, _transforms, _masks,
                                     navigation.nav_links());
                sfi.index = obj_idx;
                sfi.link = obj.edge()[0];

                // Update the intersection with a new one
                (*navigation.kernel.next) = sfi;
                navigation.distance_to_next = sfi.path;

                if (sfi.status == e_inside) {
                    // We may be on next object (trust level is high)
                    if (std::abs(sfi.path) < navigation.on_object_tolerance) {
                        navigation.object_index = obj_idx;
                        navigation.status = e_on_object;
                        navigation.trust_level = e_high_trust;
                        // Set the next object we want to reach might be end()
                        ++navigation.kernel.next;
                    }
                    // we are certainly not on the next object
                    else {
                        navigation.status = e_towards_object;
                        // Trust fully again
                        navigation.trust_level = e_full_trust;
                    }
                    // Don't sort again when coming from high trust
                    return;
                }
                // If not inside: increase and switch to next
                ++navigation.kernel.next;
            }
        }
        // Loop over all candidates and intersect again all candidates
        // - do this when your trust level is low
        else if (navigation.trust_level == e_fair_trust/* or
                 is_exhausted(navigation.kernel)*/) {
            for (auto &candidate : navigation.kernel.candidates) {
                dindex obj_idx = candidate.index;
                auto &obj = _objects[obj_idx];
                auto sfi = intersect(track, obj, _transforms, _masks,
                                     navigation.nav_links());
                candidate = sfi;
                candidate.index = obj_idx;
                candidate.link = obj.edge()[0];
            }
            set_next(navigation);
            return;
        }
        // If we end here, something went seriously wrong
        navigation.abort();
    }

    /** Helper method to sort within the kernel and find next object
     *
     * @param navigation [in, out] navigation state that contains the kernel
     */
    void set_next(state &navigation) const {

        auto &kernel = navigation.kernel;
        // Sort distance to next & set navigation status
        if (not kernel.candidates.empty()) {

            // Generally, we are on our way to some candidate
            navigation.status = e_towards_object;
            // This is only called after full (re-)evaluation
            navigation.trust_level = e_full_trust;

            // Take the nearest candidate first
            std::sort(kernel.candidates.begin(), kernel.candidates.end());
            kernel.next = kernel.candidates.begin();

            // Are we still on an object from a previous navigation pass? Then
            // goto the next candidate.
            // This also excludes adjacent portals -> we are towards obj again
            if (navigation.distance_to_next < navigation.on_object_tolerance) {
                // The object we are on
                navigation.object_index = kernel.next->index;
                // The next object that we want to reach
                ++navigation.kernel.next;
                // We might be wrong in this assumption, re-evaluate distance
                // to next in the next pass
                navigation.trust_level = e_high_trust;
            }
            // No current object
            else {
                navigation.object_index = dindex_invalid;
            }

            navigation.distance_to_next = kernel.next->path;
            return;
        }
    }

    /** Helper method to check and perform a volume switch
     *
     * @param navigation is the navigation state
     *
     * @return a flag if the volume navigation still has a heartbeat
     */
    void check_volume_switch(state &navigation) const {
        // Check if we need to switch volume index and (re-)initialize
        if (navigation.status == e_on_object and
            navigation.volume_index != navigation.current()->link) {

            // Set volume index to the next volume provided by the object
            navigation.volume_index = navigation.current()->link;
            navigation.kernel.clear();
            navigation.trust_level = e_no_trust;

            // We reached the end of the detector world
            if (navigation.volume_index == dindex_invalid) {
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
    /// the containers for all data
    const volume_container &_volumes;
    const object_container &_objects;
    const transform_container &_transforms;
    const mask_container &_masks;
};

}  // namespace detray