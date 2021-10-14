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
class single_state_navigator {

    public:
    using object = typename object_container::value_type;
    using link = typename object::edge_links;

    /** Navigation status flag */
    enum navigation_status : int {
        // e_on_target = -3,
        e_abort = -2,
        e_unknown = -1,
        e_towards_target = 0,
        e_on_target = 1,
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
        const object *on = nullptr;

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
        friend class single_state_navigator;

        public:
        /** Scalar representation of the navigation state,
         * @returns distance to next
         **/
        scalar operator()() const { return distance_to_next; }

        /** Current candidates */
        const auto &nav_kernel() { return kernel; }

        /** Current volume */
        const auto &volume() { return volume_index; }

        /** Tolerance to determine of we are on target */
        const auto &tolerance() { return on_target_tolerance; }

        /** Adjust the on-target tolerance */
        const auto &set_tolerance(scalar tol) { on_target_tolerance = tol; }

        /** get the navigation inspector */
        const auto &nav_inspector() { return inspector; }

        /** Current navigation status */
        const auto &nav_status() { return status; }

        /** Current object the navigator is on (might be invalid if between
         * objects)
         */
        const auto &on_object() { return target_index; }

        /** The links (next volume, next target finder) of current
         * candidate
         */
        const auto &nav_links() { return links; }

        /** Navigation trust level */
        const auto &nav_trust_level() { return trust_level; }

        private:
        /** Navigation state cannot be recovered from. Leave the other
         * data for inspection.
         *
         * @return navigation heartbeat
         */
        bool abort() {
            status = e_abort;
            trust_level = e_no_trust;
            return false;
        }

        /** Kernel for the targets */
        navigation_kernel<> kernel;

        /** Volume we are currently navigating in */
        dindex volume_index = dindex_invalid;

        /**  Distance to next - will be cast into a scalar with call operator
         */
        scalar distance_to_next = std::numeric_limits<scalar>::infinity();

        /** The on target tolerance - permille */
        scalar on_target_tolerance = 1e-3;

        /** The inspector type of this navigation engine */
        inspector_type inspector;

        /**  The navigation status */
        navigation_status status = e_unknown;

        /** Index of a target (surface/portal) if is reached, otherwise
         * invalid
         */
        dindex target_index = dindex_invalid;

        // Point to the next volume and target finder
        link links = {};

        /** The navigation trust level */
        navigation_trust_level trust_level = e_no_trust;
    };

    /** Constructor with move constructor
     *
     * @param data the container for all data: volumes, primitives (objects
     *  i.e. targets and portals), transforms and masks (in a tuple by
     *  type)
     */
    single_state_navigator(const volume_container &volumes,
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
    bool status(state &navigation, const track_t &track) const {

        bool heartbeat = true;

        // Retrieve the volume & set index.
        const auto &volume = _volumes[navigation.volume_index];

        // If there is no_trust (e.g. at the beginning of navigation), the
        // kernel will be initialized. Otherwise the candidates are
        // re-evaluated based on current trust level
        update_kernel(navigation, track, volume.range());
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
        bool heartbeat = true;

        // FWe are already on the right track, nothing left to do
        if (navigation.trust_level == e_full_trust) {
            return heartbeat;
        }

        // Retrieve the current volume
        const auto &volume = _volumes[navigation.volume_index];

        if (navigation.kernel.empty()) {
            return navigation.abort();
        }
        if (is_exhausted(navigation.kernel)) {
            navigation.kernel.clear();
            navigation.trust_level = e_no_trust;
        }
        update_kernel(navigation, track, volume.range());
        navigation.inspector(navigation);
        return heartbeat;
    }

    /** Helper method to intersect all objects of a target/portal store
     *
     * @tparam range the type of range in the detector data containers
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param track the track information
     * @param obj_range the target/portal index range in the detector cont
     * @param on_object ignores on target solution
     *
     */
    template <typename track_t, typename range_t>
    void initialize_kernel(state &navigation, const track_t &track,
                           const range_t &obj_range,
                           bool on_object = false) const {

        // Get the number of candidates & run them through the kernel
        navigation.kernel.candidates.reserve(obj_range[1] - obj_range[0]);

        // Loop over all indexed targets, intersect and fill
        // @todo - will come from the local object finder
        for (size_t obj_idx = obj_range[0]; obj_idx < obj_range[1]; obj_idx++) {
            // Get next target
            const auto &object = _objects[obj_idx];
            // Retrieve candidate from the target
            auto sfi = intersect(track, object, _transforms, _masks,
                                 navigation.links());
            // Candidate is invalid if it oversteps too far (this is neg!)
            if (sfi.path < track.overstep_tolerance) {
                continue;
            }
            // Accept if inside, but not if the same object is excluded
            if (sfi.status == e_inside and
                (/*not on_object or*/
                 std::abs(sfi.path) > navigation.on_target_tolerance)) {
                navigation.status = e_towards_target;

                // target the candidate belongs to
                sfi.index = obj_idx;
                // the next volume if we encounter the candidate
                sfi.link = navigation.links()[0];
                navigation.kernel.candidates.push_back(sfi);
            }
        }
        // Prepare for evaluation of candidates
        sort_and_set(navigation, navigation.kernel);
    }

    /** Helper method to the update the next candidate intersection
     *
     * @tparam range the type of range in the detector data containers
     *
     * @param navigation [in, out] navigation state that contains the kernel
     * @param track the track information
     * @param obj_range the target/portal index range in the detector cont
     *
     * @return A boolean condition if kernel is exhausted or not
     */
    template <typename track_t, typename range_t>
    void update_kernel(state &navigation, const track_t &track,
                       const range_t &obj_range) const {

        // Update current candidate, or step further
        // - do this only when you trust level is high
        if (navigation.trust_level >= e_high_trust) {
            while (not is_exhausted(navigation.kernel)) {
                // Only update the last intersection
                dindex obj_idx = navigation.kernel.next->index;
                const auto &obj = _objects[obj_idx];
                auto sfi = intersect(track, obj, _transforms, _masks,
                                     navigation.links());
                sfi.index = obj_idx;
                sfi.link = navigation.links()[0];
                if (sfi.status == e_inside) {
                    // Update the intersection with a new one
                    (*navigation.kernel.next) = sfi;
                    navigation.distance_to_next = sfi.path;

                    // We may be on target (trust level is high)
                    if (std::abs(sfi.path) < navigation.on_target_tolerance) {
                        navigation.target_index = obj_idx;
                        //++navigation.kernel.next;
                        navigation.status = e_on_target;
                        navigation.trust_level = e_high_trust;

                        // Did we hit a portal? (kernel needs to be
                        // re-initialized next time)
                        check_volume_switch(navigation);
                    }
                    // we are certainly not on target
                    else {
                        navigation.status = e_towards_target;
                        // Trust fully again
                        navigation.trust_level = e_full_trust;
                    }
                    return;
                }
                // If not successful: increase and switch to next
                ++navigation.kernel.next;
            }
        }
        // Loop over all candidates and intersect again all candidates
        // - do this when your trust level is low
        else if (navigation.trust_level == e_fair_trust or
                 is_exhausted(navigation.kernel)) {
            for (auto &candidate : navigation.kernel.candidates) {
                dindex obj_idx = candidate.index;
                auto &obj = _objects[obj_idx];
                auto sfi = intersect(track, obj, _transforms, _masks,
                                     navigation.links());
                candidate = sfi;
                candidate.index = obj_idx;
                candidate.link = navigation.links()[0];
            }
            sort_and_set(navigation);

            return;
        }
        // This kernel cannot be trusted
        navigation.trust_level = e_no_trust;
        initialize_kernel(navigation, track, obj_range);
    }

    /** Helper method to sort within the kernel
     *
     * @param navigation [in, out] navigation state that contains the kernel
     */
    void sort_and_set(state &navigation) const {

        auto kernel = navigation.kernel;
        // Sort and set distance to next & navigation status
        if (not kernel.candidates.empty()) {
            navigation.trust_level = e_full_trust;
            std::sort(kernel.candidates.begin(), kernel.candidates.end());

            // beginning of navigation with this kernel
            kernel.next = kernel.candidates.begin();
            navigation.distance_to_next = kernel.next->path;

            // Are we already/still on target? Then target the next cand.
            if (navigation.distance_to_next < navigation.on_target_tolerance) {
                navigation.status = e_on_target;
                navigation.target_index = kernel.next->index;
                check_volume_switch(navigation);
            }
            // No current object
            else {
                navigation.status = e_towards_target;
                navigation.target_index = dindex_invalid;
            }

            return;
        }

        // If after full evaluation no candidates are there, abort
        navigation.abort();
    }

    /** Helper method to check and perform a volume switch
     *
     * @param navigation is the navigation state
     *
     * @return a flag if the volume navigation still has a heartbeat
     */
    void check_volume_switch(state &navigation) const {
        // Check if we need to switch volume index and (re-)initialize
        if (navigation.status == e_on_target and
            navigation.volume_index != navigation.kernel.next->link) {
            // Set volume index to the next volume provided by the target
            navigation.volume_index = navigation.kernel.next->link;
            navigation.kernel.clear();
            navigation.trust_level = e_no_trust;
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