/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <map>
#include <sstream>
#include <string>

#include "core/transform_store.hpp"
#include "geometry/index_geometry.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/populator.hpp"
#include "grids/serializer2.hpp"
#include "tools/local_object_finder.hpp"

namespace detray {

// Algebra, point2 is not strongly typed
using point3 = __plugin::point3;
using vector3 = __plugin::vector3;
using point2 = __plugin::point2;

/** The detector definition.
 *
 * This class is a heavy templated detector definition class, that sets the
 * interface between geometry, navigator and grid.
 *
 * @tparam array_type the type of the internal array, must have STL semantics
 * @tparam tuple_type the type of the internal tuple, must have STL semantics
 * @tparam vector_type the type of the internal array, must have STL semantics
 * @tparam alignable_store the type of the transform store
 * @tparam geometry_type the geometry implementation to be used
 * @tparam surfaces_populator_type the type of populator used to fill the
 * surfaces grids
 * @tparam surfaces_serializer_type the type of the memory serializer for the
 * surfaces grids
 *
 */
template <template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple,
          template <typename> class vector_type = dvector,
          typename alignable_store = static_transform_store<vector_type>,
          typename geometry_type = index_geometry<array_type, vector_type,
                                                  tuple_type, dindex, dindex>,
          typename surfaces_populator_type =
              attach_populator<false, dindex, vector_type>,
          typename surfaces_serializer_type = serializer2,
          typename name_map = std::map<dindex, std::string>>
class detector {

    public:
    /// Forward the alignable container and context
    using transform_store = alignable_store;
    using context = typename alignable_store::context;

    /// Export geometry types

    /// Forward the alignable container and context
    // using transform_store = typename geometry_type::transform_store;
    // using context = typename geometry_type::transform_store::context;

    /// Set the geometry types
    using geometry = geometry_type;
    // enumeration of the geometry objects in the geometry (surfaces/portals
    // etc.)
    using objects = typename geometry_type::known_objects;
    using mask_id = typename geometry_type::mask_id;
    // geometry oject types
    using volume = typename geometry_type::volume_type;
    using portal = typename geometry_type::portal;
    using surface = typename geometry_type::surface;

    // Determined by the geometry, due to potentially different linking in masks
    using mask_container = typename geometry_type::mask_container;

    // Source links of geometry objects, if detray acts as plugin
    using surface_source_link = typename geometry_type::source_link;
    using bounds_source_link = typename geometry_type::bounds_link;

    /// Accelerator structure

    /// Volume grid definition
    using volume_grid =
        grid2<replace_populator<dindex, std::numeric_limits<dindex>::max(),
                                vector_type>,
              axis::irregular<array_type, vector_type>,
              axis::irregular<array_type, vector_type>, serializer2>;

    using surfaces_regular_axis = axis::regular<array_type>;
    using surfaces_circular_axis = axis::circular<array_type>;
    using surfaces_regular_circular_grid =
        grid2<surfaces_populator_type, surfaces_regular_axis,
              surfaces_circular_axis, surfaces_serializer_type, array_type,
              tuple_type, vector_type>;

    // Neighborhood finder, using accelerator data structure
    using surfaces_finder = surfaces_regular_circular_grid;

    // Temporary container struct, used to fill the detector data
    using transform_container =
        array_type<transform_store, mask_id::e_mask_types>;

    detector() = delete;

    /** Allowed costructor
     * @param name the detector name
     */
    detector(const std::string &name) : _name(_name) {}

    /** Add a new volume and retrieve a reference to it
     *
     * @param bounds of the volume, they are expected to be already attaching
     * @param surfaces_finder_entry of the volume, where to entry the surface
     * finder
     *
     * @return non-const reference of the new volume
     */
    volume &new_volume(const array_type<scalar, 6> &bounds,
                       dindex surfaces_finder_entry = dindex_invalid) {

        return _geometry.new_volume(bounds, surfaces_finder_entry);
    }

    /** @return the name of the detector */
    const std::string &name() const { return _name; }

    /** @return the contained volumes of the detector - const access */
    decltype(auto) volumes() const { return _geometry.volumes(); }

    /** @return the volume by @param volume_index - non-const access */
    decltype(auto) volume_by_index(dindex volume_index) {
        return _geometry.volume_by_index(volume_index);
    }

    /** @return the volume by @param volume_index - const access */
    decltype(auto) volume_by_index(dindex volume_index) const {
        return _geometry.volume_by_index(volume_index);
    }

    /** @return the volume by @param position - const access */
    decltype(auto) volume_by_pos(const point3 &p) const {
        point2 p2 = {getter::perp(p), p[2]};
        dindex volume_index = _volume_grid.bin(p2);
        return _geometry.volume_by_index(volume_index);
    }

    /** @return all objects of a given type */
    template <objects object_type>
    decltype(auto) get_objects() const {
        return _geometry.template objects<object_type>();
    }

    /** Convenience function:
     * @return all surfaces in the detector */
    decltype(auto) surfaces() const {
        return _geometry.template objects<objects::e_surface>();
    }

    /** Convenience function:
     * @return all portals in the detector */
    decltype(auto) portals() const {
        return _geometry.template objects<objects::e_portal>();
    }

    /** Get @return all surface/portal masks in the geometry */
    decltype(auto) masks() const { return _masks; }

    /** Get the current number of transforms in the store
     *
     * @param ctx The context of the call
     *
     * @return transform store size for a given context
     */
    const unsigned int n_transforms(const context &ctx) const {
        return _transforms.size(ctx);
    }

    /** Get all transform in an index range from the detector
     *
     * @param range The range of surfaces/portals in the transform store
     * @param ctx The context of the call
     *
     * @return ranged iterator to the object transforms
     */
    const auto transforms(const dindex_range range,
                          const context &ctx = {}) const {
        return _transforms.range(std::move(range), ctx);
    }

    /** Get all transform in an index range from the detector
     *
     * @param range The range of surfaces/portals in the transform store
     * @param ctx The context of the call
     *
     * @return ranged iterator to the object transforms
     */
    const auto &transforms(const context &ctx = {}) const {
        return _transforms;
    }

    /** Add a new full set of alignable transforms and masks for geometry
     * objects. These are added to the detector in one go to ensure constitent
     * indexing.
     *
     * @tparam current_type the current mask context to be processed
     * @tparam object_container surfaces/portals for which the links are updated
     * @tparam mask_container surface/portal masks, sorted by type
     * @tparam object_type check whether we deal with surfaces or portals
     *
     * @param volume The volume we add the transforms to
     * @param objects The geometry objects in the volume
     * @param trfs The transforms
     * @param ctx The context of the call
     *
     * @note can throw an exception if input data is inconsistent
     */
    template <objects object_type, typename object_container,
              typename mask_container, typename transform_container>
    inline void add_objects(
        volume &volume, object_container &objects, mask_container &masks,
        transform_container &trfs,
        const typename alignable_store::context ctx = {}) noexcept(false) {
        unroll_container_filling<0, object_container, mask_container,
                                 transform_container, object_type>(
            volume, objects, masks, trfs, ctx);
    }

    /** Convenience function to add detector surfaces. See
     * detector::add_objects(). */
    template <typename object_container>
    inline void add_surfaces(
        volume &volume, object_container &surfaces, mask_container &masks,
        transform_container &trfs,
        const typename alignable_store::context ctx = {}) noexcept(false) {
        add_objects<objects::e_surface>(volume, surfaces, masks, trfs, ctx);
    }

    /** Convenience function to add detector portals. See
     * detector::add_objects(). */
    template <typename object_container>
    inline void add_portals(
        volume &volume, object_container &portals, mask_container &masks,
        transform_container &trfs,
        const typename alignable_store::context ctx = {}) noexcept(false) {
        add_objects<objects::e_portal>(volume, portals, masks, trfs, ctx);
    }

    /** Unrolls the data containers according to the mask type and fill the
     *  global containers. Register the indexing in the geometry.
     *
     * @tparam current_type the current mask context to be processed
     * @tparam object_container surfaces/portals for which the links are updated
     * @tparam mask_container surface/portal masks, sorted by type
     * @tparam object_type check whether we deal with surfaces or portals
     *
     * @param volume The volume we add the transforms to
     * @param objects The geometry objects in the volume
     * @param trfs The transforms
     * @param ctx The context of the call
     *
     * @note can throw an exception if input data is inconsistent
     */
    template <size_t current_type = 0, typename object_container,
              typename mask_container, typename transform_container,
              objects object_type>
    inline void unroll_container_filling(
        volume &volume, object_container &objects, mask_container &masks,
        transform_container &trfs,
        const typename alignable_store::context ctx = {}) noexcept(false) {
        // Get the surfaces/portals for a mask type
        auto &typed_objects = objects[current_type];
        // Get the corresponding transforms
        const auto &object_transforms = trfs[current_type];
        // and the corresponding masks
        auto &object_masks = masks.template group<current_type>();

        if (object_transforms.size(ctx) != 0 and not typed_objects.empty()) {
            // Current offsets into detectors containers
            const auto trsf_offset = n_transforms(ctx);
            const auto mask_offset = _masks.template size<current_type>();

            // Fill the correct mask type
            _masks.template add_masks<current_type>(object_masks);
            _transforms.append(ctx, std::move(std::get<current_type>(trfs)));

            // Update the surfaces mask link
            for (auto &obj : typed_objects) {
                _geometry.update_mask_link(obj, mask_offset);
                _geometry.update_transform_link(obj, trsf_offset);
            }

            // Now put the updates objects into the geometry
            _geometry.template add_objects<object_type>(volume, typed_objects);
        }

        // Next mask type
        if constexpr (current_type <
                      std::tuple_size_v<typename mask_container::mask_tuple> -
                          1) {
            return unroll_container_filling<current_type + 1, object_container,
                                            mask_container, transform_container,
                                            object_type>(volume, objects, masks,
                                                         trfs, ctx);
        }
        // If no mask type fits, don't fill the data.
    }

    /** Add the volume grid - move semantics
     *
     * @param v_grid the volume grid to be added
     */
    void add_volume_grid(volume_grid &&v_grid) {
        _volume_grid = std::move(v_grid);
    }

    /** @return the volume grid - const access */
    const volume_grid &volume_search_grid() const { return _volume_grid; }

    /** Add local surface finders linked to from the portals - move semantics
     *
     * This connects portals and surface grids
     */
    void add_surfaces_finders(vector_type<surfaces_finder> &&surfaces_finders) {
        _surfaces_finders = std::move(surfaces_finders);
    }

    /** @return the surface finders - const access */
    const vector_type<surfaces_finder> &surfaces_finders() const {
        return _surfaces_finders;
    }

    /** Output to string */
    const std::string to_string(const name_map &names) const {
        std::stringstream ss;

        ss << "[>] Detector '" << _name << "' has " << _geometry.n_volumes()
           << " volumes." << std::endl;
        ss << "    contains  " << _surfaces_finders.size()
           << " local surface finders." << std::endl;

        for (const auto &[i, v] : enumerate(_geometry.volumes())) {
            ss << "[>>] Volume at index " << i << ": " << std::endl;
            ss << " - name: '" << v.name(names) << "'" << std::endl;

            ss << "     contains    " << v.template n_objects<true>() << " surfaces "
               << std::endl;

            ss << "                 " << v.template n_objects<false>() << " portals "
               << std::endl;

            if (v.surfaces_finder_entry() != dindex_invalid) {
                ss << "  sf finders idx " << v.surfaces_finder_entry() << std::endl;
            }

            const auto& bounds = v.bounds();
            ss << "     bounds r = (" << bounds[0] << ", " << bounds[1] << ")"
               << std::endl;
            ss << "            z = (" << bounds[2] << ", " << bounds[3] << ")"
               << std::endl;
        }

        return ss.str();
    };

    private:
    std::string _name = "unknown_detector";

    /** Keeps the geometry object linking*/
    geometry_type _geometry = {};

    /** Keeps all of the transform data in contiguous memory*/
    transform_store _transforms = {};

    /** Surface and portal masks of the detector in contigous memory */
    mask_container _masks = {};

    vector_type<surfaces_finder> _surfaces_finders;

    volume_grid _volume_grid = volume_grid(std::move(axis::irregular{{}}),
                                           std::move(axis::irregular{{}}));
};

}  // namespace detray
