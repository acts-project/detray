/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <sstream>
#include <string>

#include "core/intersection.hpp"
#include "core/surface_base.hpp"
#include "core/transform_store.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/populator.hpp"
#include "grids/serializer2.hpp"
#include "masks/masks.hpp"
#include "tools/concentric_cylinder_intersector.hpp"
#include "tools/cylinder_intersector.hpp"
#include "tools/local_object_finder.hpp"
#include "tools/planar_intersector.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {

// Algebra, point2 is not strongly typed
using point3 = __plugin::point3;
using vector3 = __plugin::vector3;
using point2 = __plugin::point2;

/** Indexed detector definition.
 *
 * This class is a heavy templated detector definition class, that connects
 * surfaces, layers and volumes via an indexed system.
 *
 * @tparam array_type the type of the internal array, must have STL semantics
 * @tparam tuple_type the type of the internal tuple, must have STL semantics
 * @tparam vector_type the type of the internal array, must have STL semantics
 * @tparam alignable_store the type of the transform store
 * @tparam surface_source_link the type of the link to an external surface
 * source
 * @tparam bounds_source_link the type of the link to an external bounds source
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
          typename surface_source_link = dindex,
          typename bounds_source_link = dindex,
          typename surfaces_populator_type =
              attach_populator<false, dindex, vector_type>,
          typename surfaces_serializer_type = serializer2>
class detector {

    public:
    /** Encodes the position in a collection container for the respective
        mask type (not for portals for the moment). */
    enum mask_container_index : unsigned int {
        e_mask_types = 5,
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_ring2 = 4,
        e_single3 = std::numeric_limits<unsigned int>::max(),
        e_unknown = std::numeric_limits<unsigned int>::max(),
    };

    enum use_primitive : bool {
        e_surface = true,
        e_portal = false,
    };

    /// The detector type
    using detector_type = detector<array_type, tuple_type, vector_type, alignable_store, surface_source_link, bounds_source_link, surfaces_populator_type, surfaces_serializer_type>;

    /// Forward the alignable container and context
    using transform_store = alignable_store;
    using context = typename alignable_store::context;

    /// Volume grid definition
    using volume_grid = grid2<replace_populator<dindex, std::numeric_limits<dindex>::max(), vector_type>,
                                  axis::irregular<array_type, vector_type>,
                                  axis::irregular<array_type, vector_type>,
                                  serializer2>;

    /// Portals components:
    /// - links:  next volume, next (local) object finder
    using volume_links = array_type<dindex, 2>;

    /// - mask types
    using rectangle = rectangle2<planar_intersector, __plugin::cartesian2, volume_links, 0>;
    using trapezoid = trapezoid2<planar_intersector, __plugin::cartesian2, volume_links, 1>;
    using annulus = annulus2<planar_intersector, __plugin::cartesian2, volume_links, 2>;
    using cylinder = cylinder3<false, cylinder_intersector, __plugin::cylindrical2, volume_links, 3>;
    using disc = ring2<planar_intersector, __plugin::cartesian2, volume_links, 4>;

    // - mask index: type, { first/last }
    using portal_mask_index = tuple_type<dindex, array_type<dindex, 2>>;

    /** The Portal definition:
     *  <transform_link, mask_index, volume_link, source_link >
     * 
     * transform_link: index into the transform container
     * mask_index: typed index into the mask container
     * volume_link: index of the volume this portal belongs to
     * source_link: some link to an eventual exernal representation
     * 
     */
    using portal = surface_base<dindex, portal_mask_index, dindex, surface_source_link>;
    using portal_container = vector_type<portal>;

    /// Surface components:
    /// - masks, with mask identifiers 0,1,2

    /// - mask index: type, entry
    using surface_mask_index = array_type<dindex, 2>;
    using mask_container = tuple_type<vector_type<rectangle>,
                                      vector_type<trapezoid>,
                                      vector_type<annulus>,
                                      vector_type<cylinder>,
                                      vector_type<disc>>;

    using surface_link = surface_source_link;
    /** The Surface definition:
     *  <transform_link, mask_link, volume_link, source_link >
     */
    using surface = surface_base<dindex, surface_mask_index, dindex, surface_link>;
    using surface_container = vector_type<surface>;

    using surfaces_regular_axis = axis::regular<array_type>;
    using surfaces_circular_axis = axis::circular<array_type>;
    using surfaces_regular_circular_grid = grid2<surfaces_populator_type,
                                                 surfaces_regular_axis,
                                                 surfaces_circular_axis,
                                                 surfaces_serializer_type,
                                                 array_type,
                                                 tuple_type,
                                                 vector_type>;

    using surfaces_finder = surfaces_regular_circular_grid;

    /** Temporary container structures that are used to fill the detector. 
     * The respective objects are sorted by mask type, so that they can be 
     * unrolled and filled in lockstep with the masks
     */
    using surface_filling_container = array_type<vector_type<surface>, 5>;
    using portal_filling_container = array_type<vector_type<portal>, 5>;
    using transform_container = array_type<transform_store, 5>;

    /** Nested volume struct that holds the local information of the
     * volume and its portals.
     */
    class volume {
        friend class detector<array_type, tuple_type, vector_type,
                              alignable_store, surface_source_link,
                              bounds_source_link, surfaces_populator_type,
                              surfaces_serializer_type>;

        public:
        /** Deleted constructor */
        volume() = delete;


        /** The Surface definition:
         *  <transform_link, mask_link, volume_link, source_link >
         */
        volume(const std::string &name) : _name(name) {}

        /** Contructor with name and bounds
         * @param name of the volume
         * @param bounds of the volume
         * @param d detector the volume belongs to
         */
        volume(const std::string &name, const array_type<scalar, 6> &bounds)
            : _name(name), _bounds(bounds) {}

        /** Copy ctor makes sure constituents keep valid volume pointer
         *
         * @param other Volume to be copied
         */
        volume(const volume &other) = default;

        /** @return the bounds - const access */
        const array_type<scalar, 6> &bounds() const { return _bounds; }

        /** @return the name */
        const std::string &name() const { return _name; }

        /** @return the index */
        dindex index() const { return _index; }

        /** @return the entry into the local surface finders */
        dindex surfaces_finder_entry() const { return _surfaces_finder_entry; }

        /** @return if the volume is empty or not */
        bool empty() const { return is_empty_range(_surface_range); }

        /** @return the number of surfaces in the volume */
        template <bool use_surfaces = e_surface>
        dindex n_objects() {
            if constexpr (use_surfaces) {
                return n_in_range(_surface_range);
            } else {
                return n_in_range(_portal_range);
            }
        }

        /** @return the number of surfaces in the volume */
        template <bool use_surfaces = e_surface>
        const dindex n_objects() const {
            if constexpr (use_surfaces) {
                return n_in_range(_surface_range);
            } else {
                return n_in_range(_portal_range);
            }
        }

        /** Set the index into the detector surface container
         *
         * @param range Surface index range
         */
        template <bool surface_range = e_surface>
        void set_range(dindex_range range) {
            if constexpr (surface_range) {
                update_range(_surface_range, std::move(range));
            } else {
                update_range(_portal_range, std::move(range));
            }
        }

        /** @return range of surfaces- const access */
        template <bool surface_range = e_surface>
        const auto &range() const {

            if constexpr (surface_range) {
                return _surface_range;
            } else {
                return _portal_range;
            }
        }

        /** @return range of portals - const access */
        // const auto &portal_range() const { return _portal_range; }

        /** @return range of surface transforms - const access */
        template <bool surface_range = e_surface>
        const auto &trf_range() const {
            if constexpr (surface_range) {
                return _surface_trf_range;
            } else {
                return _portal_trf_range;
            }
        }

        /** Set the index into the detector transform store for portals
         *
         * @param range Portal transform index range
         */
        template <bool surface_range = e_surface>
        void set_trf_range(dindex_range range) {
            if constexpr (surface_range) {
                update_range(_surface_trf_range, std::move(range));
            } else {
                update_range(_portal_trf_range, std::move(range));
            }
        }

        private:
        /** Volume section: name */
        std::string _name = "unknown";

        /** Volume index */
        dindex _index = dindex_invalid;

        /** Transform ranges in the detector transform store.*/
        dindex_range _surface_trf_range = {dindex_invalid, dindex_invalid};
        dindex_range _portal_trf_range = {dindex_invalid, dindex_invalid};

        /** Index ranges in the detector surface/portal containers.*/
        dindex_range _surface_range = {dindex_invalid, dindex_invalid};
        dindex_range _portal_range = {dindex_invalid, dindex_invalid};

        /** Index into the surface finder container */
        dindex _surfaces_finder_entry = dindex_invalid;

        /** Bounds section, default for r, z, phi */
        array_type<scalar, 6> _bounds = {0.,
                                         std::numeric_limits<scalar>::max(),
                                         -std::numeric_limits<scalar>::max(),
                                         std::numeric_limits<scalar>::max(),
                                         -M_PI,
                                         M_PI};

        /**
         * @param range Any index range
         *
         * @return the number of indexed objects
         */
        inline dindex n_in_range(const dindex_range &range) {
            return range[1] - range[0];
        }

        /**
         * @param range Any index range
         *
         * @return the number of indexed objects
         */
        inline const dindex n_in_range(const dindex_range &range) const {
            return range[1] - range[0];

        /** Test whether a range is empty
         *
         * @param range Any index range
         *
         * @return boolean whether the range is empty
         */
        inline bool is_empty_range(const dindex_range &range) {
            return n_in_range(range) == 0;
        }

        /** Test whether a range is empty - const
         *
         * @param range Any index range
         *
         * @return boolean whether the range is empty
         */
        inline const bool is_empty_range(const dindex_range &range) const {
            return n_in_range(range) == 0;
        }

        /** Set or update a range
         *
         * @param range One of the volume member ranges
         * @param other new index range
         *
         * @return boolean whether the range is empty
         */
        inline void update_range(dindex_range &range, dindex_range &&other) {
            // Range not set yet
            if (range[0] == dindex_invalid) {
                range = other;
            } else {
                range[1] += other[1] - other[0];
            }
        }
    };

    /** Allowed costructor
     * @param name the detector
     */
    detector(const std::string &name) : _name(name) {}

    /** Copy constructor makes sure the volumes belong to new detector.
     *
     * @param other Detector to be copied
     */
    detector(const detector &other) = default;
    detector() = delete;
    ~detector() = default;

    /** Add a new volume and retrieve a reference to it
     *
     * @param name of the volume
     * @param bounds of the volume, they are expected to be already attaching
     * @param surfaces_finder_entry of the volume, where to entry the surface
     * finder
     *
     * @return non-const reference of the new volume
     */
    volume &new_volume(const std::string &name,
                       const array_type<scalar, 6> &bounds,
                       dindex surfaces_finder_entry = dindex_invalid) {
        _volumes.emplace_back(name, bounds);
        dindex cvolume_idx = _volumes.size() - 1;
        volume &cvolume = _volumes[cvolume_idx];
        cvolume._index = cvolume_idx;
        cvolume._surfaces_finder_entry = surfaces_finder_entry;
        return cvolume;
    }

    /** @return the name of the detector */
    const std::string &name() const { return _name; }

    /** @return the contained volumes of the detector - const access */
    const vector_type<volume> &volumes() const { return _volumes; }

    /** @return the volume by @param volume_index - const access */
    const volume &indexed_volume(dindex volume_index) const {
        return _volumes[volume_index];
    }

    /** @return the volume by @param position - const access */
    const volume &indexed_volume(const point3 &p) const {
        point2 p2 = {getter::perp(p), p[2]};
        dindex volume_index = _volume_grid.bin(p2);
        return _volumes[volume_index];
    }

    /** @return the volume by @param volume_index - non-const access */
    volume &indexed_volume(dindex volume_index) {
        return _volumes[volume_index];
    }

    /** Add a new full set of alignable transforms for surfaces of a volume
     *
     * @param volume The volume we add the transforms to
     * @param surfaces The surfaces in the volume
     * @param trfs The surface transforms (same number as surfaces)
     * @param ctx The context of the call
     *
     * @note can throw an exception if input data is inconsistent
     */
    template <bool add_surfaces = e_surface, typename object_container,
              typename mask_container>
    void add_objects(
        volume &volume, object_container &surfaces, mask_container &masks,
        transform_container &trfs,
        const typename alignable_store::context ctx = {}) noexcept(false) {
        unroll_container_filling<0, object_container, mask_container,
                                 add_surfaces>(surfaces, masks, volume, trfs,
                                               ctx);
    }

    /** @return all surfaces/portals in the detector */
    template <bool get_surface = true>
    const auto &surfaces() const {
        if constexpr (get_surface) {
            return _surfaces;
        } else {
            return _portals;
        }
    }

    /** Get @return all surface/portal masks in the detector */
    template <bool surface_masks = e_surface>
    const auto &masks() const {
        if constexpr (surface_masks) {
            return _surface_masks;
        } else {
            return _portal_masks;
        }
    }

    /** Get @return all surface/portal masks in the detector */
        
        template<bool surface_masks = e_surface>
        const auto& masks() const
        {
            return _masks;
        }

    /** Unrolls the data containers according to mask type and fills the
     * detector.
     *
     * @tparam current_idx the current mask context to be processed
     * @tparam object_container surfaces/portals for which the links are updated
     * @tparam mask_container surface/portal masks, sorted by type
     * @tparam add_surfaces check whether we deal with surfaces or portals
     *
     */
    template <size_t current_idx = 0,
              typename object_container,
              typename mask_container,
              bool add_surfaces = e_surface>
    void unroll_container_filling(object_container &objects,
                                  mask_container &masks,
                                  volume &volume,
                                  transform_container &trfs,
                                  const typename alignable_store::context ctx = {})
    {
        // Get the surfaces/portals for a mask type
        auto &typed_objects = objects[current_idx];
        // Get the corresponding transforms
        const auto &surface_transforms = trfs[current_idx];
        // and the corresponding masks
        auto &object_masks = std::get<current_idx>(masks);

        // Fill object masks into the correct detector container
        auto add_detector_masks = [&](auto &detector_container) -> dindex
        {
            auto &detector_masks = std::get<current_idx>(detector_container);
            const dindex mask_offset = detector_masks.size();
            detector_masks.reserve(mask_offset + object_masks.size());
            detector_masks.insert(detector_masks.end(), object_masks.begin(), object_masks.end());

            return mask_offset;
        };

        // Fill objects (surfaces or portals) into the detector container
        auto add_detector_objects = [&](auto &detector_container)
        {
            detector_container.reserve(detector_container.size() + typed_objects.size());
            detector_container.insert(detector_container.end(), typed_objects.begin(), typed_objects.end());
        };


        if (surface_transforms.size(ctx) != 0 and not typed_objects.empty())
        {
            // Current offsets into detectors containers
            const auto trsf_offset = transform_index(ctx);
            _transforms.append(ctx, std::move(std::get<current_idx>(trfs)));

            // Fill surface or portal container?
            if constexpr (add_surfaces)
            {
                // Surface transforms for this volume
                volume.template set_trf_range<e_surface>(dindex_range{trsf_offset, transform_index(ctx)});

                const auto sf_mask_offset = add_detector_masks(_masks);

                // Update the surfaces mask link
                for (auto &sf : typed_objects)
                {
                    std::get<1>(sf.mask()) += sf_mask_offset;
                }

                // Now put the surfaces into the detector
                add_detector_objects(_surfaces);
                volume.template set_range<e_surface>({_surfaces.size() - typed_objects.size(), _surfaces.size()});
            }
            else
            {
                // Portal transforms for this volume
                volume.template set_trf_range<e_portal>(dindex_range{trsf_offset, transform_index(ctx)});

                // Fill the correct mask type
                const auto pt_mask_offset = add_detector_masks(_masks);

                // Update the portals mask links
                for (auto &obj : typed_objects)
                {
                    auto& portal_mask_index = std::get<1>(obj.mask());
                    portal_mask_index[0] += pt_mask_offset;
                    portal_mask_index[1] += pt_mask_offset;
                }

                // Now put the portals into the detector
                add_detector_objects(_portals);
                volume.template set_range<e_portal>({_portals.size() - typed_objects.size(), _portals.size()});
            }
        }
        // Next mask type
        if constexpr (current_idx < std::tuple_size_v<mask_container> - 1)
        {
            return unroll_container_filling<current_idx + 1, object_container, mask_container, add_surfaces>(objects, masks, volume, trfs, ctx);
        }
    }

    /** Get the current transform index
     *
     * @param ctx The context of the call
     *
     * @return Index to add new transforms at
     */
    const unsigned int transform_index(const context &ctx) const {
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
    const std::string to_string() const {
        std::stringstream ss;
        ss << "[>] Detector '" << _name << "' has " << _volumes.size()
           << " volumes." << std::endl;
        ss << "    contains  " << _surfaces_finders.size()
           << " local surface finders." << std::endl;
        for (const auto &[i, v] : enumerate(_volumes)) {
            ss << "[>>] Volume at index " << i << " - name: '" << v.name()
               << "'" << std::endl;
            ss << "     contains    " << v.template n_objects<e_surface>()
               << " detector surfaces" << std::endl;
            ss << "                 " << v.template n_objects<e_portal>()
               << " detector portals" << std::endl;
            if (v._surfaces_finder_entry != dindex_invalid) {
                ss << "     finders idx " << v._surfaces_finder_entry
                   << std::endl;
            }
            ss << "     bounds r = (" << v._bounds[0] << ", " << v._bounds[1]
               << ")" << std::endl;
            ss << "            z = (" << v._bounds[2] << ", " << v._bounds[3]
               << ")" << std::endl;
        }
        return ss.str();
    };

    private:
    std::string _name = "unknown_detector";

    /** Contains the geometrical relations*/
    vector_type<volume> _volumes = {};

    /** Keeps all of the transform data in contiguous memory*/
    transform_store _transforms = {};

    /** All surfaces and portals in the detector in contigous memory */
    surface_container _surfaces = {};
    portal_container _portals = {};

    /** Surface and portal masks of the detector in contigous memory */
    mask_container _masks = {};

    vector_type<surfaces_finder> _surfaces_finders;

    volume_grid _volume_grid = volume_grid(std::move(axis::irregular{{}}),
                                           std::move(axis::irregular{{}}));
};

}  // namespace detray
