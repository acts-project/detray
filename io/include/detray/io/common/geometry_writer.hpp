/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/io/common/detail/definitions.hpp"
#include "detray/io/common/detail/utils.hpp"
#include "detray/io/common/grid_writer.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/type_registry.hpp"

// System include(s)
#include <algorithm>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace detray {

/// @brief Abstract base class for tracking geometry writers
template <class detector_t>
class geometry_writer : public writer_interface<detector_t> {

    using base_type = writer_interface<detector_t>;

    protected:
    /// Tag the writer as "geometry"
    inline static const std::string tag = "geometry";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Serialize the header information into its payload
    static geo_header_payload write_header(const detector_t& det,
                                           const std::string_view det_name) {
        geo_header_payload header_data;

        header_data.version = detail::get_detray_version();
        header_data.detector = det_name;
        header_data.tag = tag;
        header_data.date = detail::get_current_date();

        header_data.n_volumes = det.volumes().size();
        header_data.n_surfaces = det.n_surfaces();

        return header_data;
    }

    /// Serialize a detector @param det into its io payload
    static detector_payload serialize(const detector_t& det) {
        detector_payload det_data;
        det_data.volumes.reserve((det.volumes().size()));

        for (const auto& vol : det.volumes()) {
            det_data.volumes.push_back(serialize(vol, det));
        }

        return det_data;
    }

    /// Serialize a link @param idx into its io payload
    static single_link_payload serialize(const std::size_t idx) {
        single_link_payload link_data;
        link_data.link = idx;

        return link_data;
    }

    /// Serialize a surface transform @param trf into its io payload
    static transform_payload serialize(
        const typename detector_t::transform3& trf) {
        transform_payload trf_data;

        const auto& t = trf.translation();
        const auto& x = trf.x();
        const auto& y = trf.y();
        const auto& z = trf.z();

        trf_data.tr = {t[0], t[1], t[2]};
        trf_data.rot = {x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2]};

        return trf_data;
    }

    /// Serialize a surface mask @param m into its io payload
    template <typename mask_t,
              std::enable_if_t<!std::is_same_v<typename mask_t::shape, void>,
                               bool> = true>
    static mask_payload serialize(const mask_t& m) {
        mask_payload mask_data;

        mask_data.shape = get_shape_id<typename mask_t::shape>();

        mask_data.volume_link = serialize(m.volume_link());

        mask_data.boundaries.resize(mask_t::boundaries::e_size);
        std::copy(std::cbegin(m.values()), std::cend(m.values()),
                  std::begin(mask_data.boundaries));

        return mask_data;
    }

    /// Serialize a surface material link @param m into its io payload
    template <class material_t>
    static material_link_payload serialize(const std::size_t idx) {
        using scalar_t = typename material_t::scalar_type;
        using type_id = material_link_payload::material_type;

        material_link_payload mat_data;

        // Find the correct material type index (use name for simplicity)
        if constexpr (std::is_same_v<material_t, material_slab<scalar_t>>) {
            mat_data.type = type_id::slab;
        } else if constexpr (std::is_same_v<material_t,
                                            material_rod<scalar_t>>) {
            mat_data.type = type_id::rod;
        } else {
            mat_data.type = type_id::unknown;
        }

        mat_data.index = idx;

        return mat_data;
    }

    /// Serialize a detector surface @param sf into its io payload
    static surface_payload serialize(const surface<detector_t>& sf) {
        surface_payload sf_data;

        sf_data.type = sf.id();
        sf_data.barcode = sf.barcode().value();
        sf_data.transform = serialize(sf.transform({}));
        sf_data.mask = sf.template visit_mask<get_mask_payload>();
        sf_data.material = sf.template visit_material<get_material_payload>();
        sf_data.source = serialize(sf.source());

        return sf_data;
    }

    /// Serialize a link @param idx into its io payload
    static acc_links_payload serialize(const acc_links_payload::acc_type id,
                                       const std::size_t idx) {
        acc_links_payload link_data;
        link_data.type = id;
        link_data.index = idx;

        return link_data;
    }

    /// Serialize a detector portal @param sf into its io payload
    static volume_payload serialize(
        const typename detector_t::volume_type& vol_desc,
        const detector_t& det) {
        volume_payload vol_data;

        vol_data.index = serialize(vol_desc.index());
        vol_data.transform =
            serialize(det.transform_store()[vol_desc.transform()]);
        vol_data.type = vol_desc.id();

        for (const auto& sf_desc : det.surface_lookup()) {
            if (sf_desc.volume() == vol_desc.index()) {
                vol_data.surfaces.push_back(serialize(surface{det, sf_desc}));
            }
        }

        // Only run the query, if object type is contained in volume
        const auto& link = vol_desc.full_link();
        // Initialize the std::optional
        if (link.size() > 1u) {
            vol_data.acc_links.emplace();
        }
        // Skip the first acceleration structure which exists in every volume
        // and is handled automatically during detector building
        for (unsigned int i = 1u; i < link.size(); ++i) {
            const auto& l = link[i];
            if (not l.is_invalid()) {
                const auto aclp =
                    det.surface_store().template visit<get_acc_link_payload>(l);
                vol_data.acc_links->push_back(aclp);
            }
        }

        return vol_data;
    }

    private:
    /// Retrieve @c mask_payload from mask_store element
    struct get_mask_payload {
        template <typename mask_group_t, typename index_t>
        inline auto operator()(const mask_group_t& mask_group,
                               const index_t& index) const {
            return geometry_writer<detector_t>::serialize(mask_group[index]);
        }
    };

    /// Retrieve @c material_link_payload from material_store element
    struct get_material_payload {
        template <typename material_group_t, typename index_t>
        inline auto operator()(const material_group_t&,
                               const index_t& index) const {
            return geometry_writer<detector_t>::template serialize<
                typename material_group_t::value_type>(index);
        }
    };

    /// Retrieve @c acc_links_payload from surface_tore collection
    struct get_acc_link_payload {
        template <typename acc_group_t, typename index_t>
        constexpr inline auto operator()(const acc_group_t&,
                                         const index_t& index) const {

            using accel_t = typename acc_group_t::value_type;

            if constexpr (detail::is_grid_v<accel_t>) {
                constexpr auto id{detail::get_grid_id<accel_t>()};

                return geometry_writer<detector_t>::serialize(id, index);
            } else {
                // This functor is only called for accelerator data structures
                // that are not 'brute force'
                return geometry_writer<detector_t>::serialize(
                    acc_links_payload::acc_type::unknown, index);
            }
        }
    };

    /// Infer the IO shape id from the shape type
    template <typename shape_t>
    static constexpr io::detail::mask_shape get_shape_id() {

        /// Register the mask shapes to the @c mask_shape enum in the io
        /// definitions
        using shape_registry =
            type_registry<detray::io::detail::mask_shape, annulus2D<>,
                          cuboid3D<>, cylinder2D<>, cylinder3D,
                          cylinder2D<false, cylinder_portal_intersector>,
                          rectangle2D<>, ring2D<>, trapezoid2D<>, line<true>,
                          line<false>, single3D<0>, single3D<1>, single3D<2>>;

        // Find the correct shape IO id;
        if constexpr (shape_registry::is_defined(shape_t{})) {
            return shape_registry::get_id(shape_t{});
        } else {
            return io::detail::mask_shape::unknown;
        }
    }
};

}  // namespace detray
