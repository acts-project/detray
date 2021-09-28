/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <sstream>
#include <string>
#include <utility>

#include "core/surface_base.hpp"
#include "core/mask_store.hpp"
#include "masks/masks.hpp"
#include "utils/enumerate.hpp"
#include "utils/indexing.hpp"

namespace detray {
/**
 * @brief Index geometry implementation
 *
 * This class provides a geometry that defines logic volumes which contain
 * the detector surfaces, joined together by dedicated portal surfaces. It
 * exports all types needed for navigation and strictly only keeps the
 * index data (links) that define the geometry relations.
 *
 * @tparam array_type the type of the internal array, must have STL
 *                    semantics
 * @tparam vector_type the type of the internal array, must have STL
 *                     semantics
 * @tparam surface_source_link the type of the link to an external surface
 *                             source
 *
 * @note The geometry knows nothing about coordinate systems. This is
 *       handeled by geometry access objects (e.g. the grid).
 */
template <template <typename, unsigned int> class array_type = darray,
          template <typename> class vector_type = dvector,
          template <typename...> class tuple_type = dtuple,
          typename surface_source_link = dindex,
          typename bounds_source_link = dindex>
class index_geometry {

    public:
    // Known primitives
    enum known_objects : bool {
        e_surface = true,
        e_portal = false,
    };

    /** Encodes the position in a collection container for the respective
        mask type . */
    enum mask_id : unsigned int {
        e_mask_types = 6,
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_portal_cylinder3 = 4,
        e_portal_ring2 = 5,
        e_single3 = std::numeric_limits<unsigned int>::max(),
        e_unknown = std::numeric_limits<unsigned int>::max(),
    };

    /// Portals components:
    /// - links:  next volume, next (local) object finder
    using portal_links = array_type<dindex, 2>;
    /// - masks, with mask identifiers 0, 1
    using portal_cylinder = cylinder3<false, cylinder_intersector,
                                      __plugin::cylindrical2, portal_links, e_portal_cylinder3>;
    using portal_disc =
        ring2<planar_intersector, __plugin::cartesian2, portal_links, e_portal_ring2>;
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
    using bounds_link = bounds_source_link;
    using portal = surface_base<dindex, portal_mask_index, dindex, bounds_link, portal_links>;
    using portal_container = vector_type<portal>;

    /// Surface components:
    /// - surface links
    using surface_links = array_type<dindex, 1>;
    /// - masks, with mask identifiers 0,1,2
    using surface_rectangle =
        rectangle2<planar_intersector, __plugin::cartesian2, surface_links,
                   e_rectangle2>;
    using surface_trapezoid =
        trapezoid2<planar_intersector, __plugin::cartesian2, surface_links,
                   e_trapezoid2>;
    using surface_annulus = annulus2<planar_intersector, __plugin::cartesian2,
                                     surface_links, e_annulus2>;
    using surface_cylinder =
        cylinder3<false, cylinder_intersector, __plugin::cylindrical2,
                  surface_links, e_cylinder3>;
    /// - mask index: type, entry
    using surface_mask_index = array_type<dindex, 2>;

    using mask_container = mask_store<tuple_type, vector_type, surface_rectangle, surface_trapezoid, surface_annulus, surface_cylinder, portal_cylinder, portal_disc>;

    using source_link = surface_source_link;
    /** The Surface definition:
     *  <transform_link, mask_link, volume_link, source_link >
     */
    using surface =
        surface_base<dindex, surface_mask_index, dindex, source_link, surface_links>;
    using surface_container = vector_type<surface>;

    /** Temporary container structures that are used to fill the geometry.
     * The respective objects are sorted by mask type, so that they can be
     * unrolled and filled in lockstep with the masks
     */
    using surface_filling_container = array_type<vector_type<surface>, e_mask_types>;
    using portal_filling_container = array_type<vector_type<portal>, e_mask_types>;

    /** Nested volume struct that holds the local information of the
     * volume and its portals.
     */
    class volume {

        public:
        /** Deleted constructor */
        volume() = delete;

        /** Allowed constructors
         * @param name of the volume
         * @param d detector the volume belongs to
         *
         * @note will be contructed boundless
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
        inline const array_type<scalar, 6> &bounds() const { return _bounds; }

        /** @return the name */
        inline const std::string &name() const { return _name; }

        /** @return the index */
        inline dindex index() const { return _index; }

        /** @param index the index */
        inline void set_index(const dindex index) { _index = index; }

        /** @return the entry into the local surface finders */
        inline dindex surfaces_finder_entry() const {
            return _surfaces_finder_entry;
        }

        /** @param entry the entry into the local surface finders */
        inline void set_surfaces_finder(const dindex entry) {
            _surfaces_finder_entry = entry;
        }

        /** @return if the volume is empty or not */
        inline bool empty() const { return is_empty_range(_surface_range); }

        /** @return the number of surfaces in the volume */
        template <bool primitive = e_surface>
        inline dindex n_objects() {
            if constexpr (primitive) {
                return n_in_range(_surface_range);
            } else {
                return n_in_range(_portal_range);
            }
        }

        /** @return the number of surfaces in the volume */
        template <bool primitive = e_surface>
        inline const dindex n_objects() const {
            if constexpr (primitive) {
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
        inline void set_range(dindex_range range) {
            if constexpr (surface_range) {
                update_range(_surface_range, std::move(range));
            } else {
                update_range(_portal_range, std::move(range));
            }
        }

        /** @return range of surfaces- const access */
        template <bool surface_range = e_surface>
        inline const auto &range() const {
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
        inline const auto &trf_range() const {
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
        inline void set_trf_range(dindex_range range) {
            if constexpr (surface_range) {
                update_range(_surface_trf_range, std::move(range));
            } else {
                update_range(_portal_trf_range, std::move(range));
            }
        }

        private:
        /** Volume section: name */
        std::string _name = "unknown";

        /** Bounds section, default for r, z, phi */
        array_type<scalar, 6> _bounds = {0.,
                                         std::numeric_limits<scalar>::max(),
                                         -std::numeric_limits<scalar>::max(),
                                         std::numeric_limits<scalar>::max(),
                                         -M_PI,
                                         M_PI};

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
        }

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

    /** @return total number of nodes (volumes) */
    const size_t n_volumes() const { return _volumes.size(); }

    /** @return all volumes in the geometry - const access. */
    const auto &volumes() const { return _volumes; }

    /** @return the volume by @param volume_index - const access. */
    inline const volume &volume_by_index(dindex volume_index) const {
        return _volumes[volume_index];
    }

    /** @return the volume by @param volume_index - non-const access. */
    inline volume &volume_by_index(dindex volume_index) {
        return _volumes[volume_index];
    }

    /** Add a new volume and retrieve a reference to it
     *
     * @param name of the volume
     * @param bounds of the volume, they are expected to be already attaching
     * @param surfaces_finder_entry of the volume, where to entry the surface
     * finder
     *
     * @return non-const reference of the new volume
     */
    inline volume &new_volume(const std::string &name,
                              const array_type<scalar, 6> &bounds,
                              dindex surfaces_finder_entry = dindex_invalid) {
        _volumes.emplace_back(name, bounds);
        dindex cvolume_idx = _volumes.size() - 1;
        volume &cvolume = _volumes[cvolume_idx];
        cvolume.set_index(cvolume_idx);
        cvolume.set_surfaces_finder(surfaces_finder_entry);
        return cvolume;
    }

    /** @return all surfaces/portals in the detector */
    template <bool get_surface = true>
    inline size_t n_objects() const {
        if constexpr (get_surface) {
            return _surfaces.size();
        } else {
            return _portals.size();
        }
    }

    /** @return all surfaces/portals in the detector */
    template <bool get_surface = true>
    inline const auto &objects() const {
        if constexpr (get_surface) {
            return _surfaces;
        } else {
            return _portals;
        }
    }

    inline void update_mask_link(surface &sf, const dindex mask_offset) {
        std::get<1>(sf.mask()) += mask_offset;
    }

    inline void update_mask_link(portal &pt, const dindex mask_offset) {
        auto &portal_mask_index = std::get<1>(pt.mask());
        portal_mask_index[0] += mask_offset;
        portal_mask_index[1] += mask_offset;
    }
    
    inline void update_transform_link(surface &sf, const dindex trsf_offset) {
        sf.transform() += trsf_offset;
    }

    inline void update_transform_link(portal &pt, const dindex trsf_offset) {
        pt.transform() += trsf_offset;
    }

    /** @return all surfaces/portals in the detector */
    inline void add_objects(volume &volume, const surface_container &surfaces) {
        const auto offset = _surfaces.size();
        _surfaces.reserve(_surfaces.size() + surfaces.size());
        _surfaces.insert(_surfaces.end(), surfaces.begin(), surfaces.end());

        volume.template set_range<e_surface>({offset, _surfaces.size()});
    }

    /** @return all surfaces/portals in the detector */
    inline void add_objects(volume &volume, const portal_container &portals) {
        const auto offset = _portals.size();
        _portals.reserve(_portals.size() + portals.size());
        _portals.insert(_portals.end(), portals.begin(), portals.end());

        volume.template set_range<e_portal>({offset, _portals.size()});
    }

    /**
     * Print geometry if an external name map is provided for the volumes.
     *
     * @param names  Lookup for the names by volume index.
     *
     * @returns the geometry description as a string
     */
    // TODO: remove names
    /*template <typename name_map>
    inline const std::string to_string(name_map &names) const*/
    inline const std::string to_string() const {
        std::stringstream ss;
        for (const auto &[i, v] : enumerate(_volumes)) {
            ss << "[>>] Volume at index " << i << " - name: '" << v.name()
               << "'" << std::endl;

            ss << "     contains    " << v.template n_objects<e_surface>()
               << " surfaces " << std::endl;

            ss << "                 " << v.template n_objects<e_portal>()
               << " portals " << std::endl;

            if (v.surfaces_finder_entry() != dindex_invalid) {
                ss << "  sf finders idx " << v.surfaces_finder_entry()
                   << std::endl;
            }
            const auto &bounds = v.bounds();
            ss << "     bounds r = (" << bounds[0] << ", " << bounds[1] << ")"
               << std::endl;
            ss << "            z = (" << bounds[2] << ", " << bounds[3] << ")"
               << std::endl;
        }
        return ss.str();
    };

    private:
    /** Contains the geometrical relations*/
    vector_type<volume> _volumes = {};

    /** All surfaces and portals in the geometry in contigous memory */
    surface_container _surfaces = {};
    portal_container _portals = {};
};

}  // namespace detray