/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/grid_reader.hpp"
#include "detray/io/common/homogeneous_material_reader.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/material_map_builder.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <stdexcept>
#include <string>

namespace detray {

namespace detail {

/// Retemplate material map builder to fit the grid builder schema
template <typename D, typename G, typename B, typename F>
using material_map_builder_t = material_map_builder<D, D::Dim, B, F>;

}  // namespace detail

/// @brief Abstract base class for surface grid readers
template <class detector_t,
          typename DIM = std::integral_constant<std::size_t, 2u>>
class material_map_reader
    : public detail::grid_reader<
          detector_t, material_slab<typename detector_t::surface_type>,
          detail::material_map_builder_t,
          std::integral_constant<std::size_t, 1u>, DIM> {

    using scalar_t = typename detector_t::surface_type;
    using grid_reader_t =
        detail::grid_reader<detector_t, material_slab<scalar_t>,
                            detail::material_map_builder_t,
                            std::integral_constant<std::size_t, 1u>, DIM>;
    using base_type = grid_reader_t;

    using material_reader_t = homogeneous_material_reader<detector_t>;

    protected:
    /// Tag the reader as "material_maps"
    inline static const std::string tag = "material_maps";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Deserialize the detector grids @param grids_data from their IO
    /// payload
    static void deserialize(
        detector_builder<typename detector_t::metadata, volume_builder>
            &det_builder,
        typename detector_t::name_map &,
        const detector_grids_payload<material_slab_payload> &grids_data) {

        using namespace axis;

        // Go through the full grid deserialization once volume material is
        // added, until then, giving the explicit bounds and binnings should be
        // enough
        using regular_t = regular<host_container_types, scalar_t>;
        using bounds2D_ts = types::list<closed<static_cast<label>(0u)>,
                                        closed<static_cast<label>(1u)>>;
        using binning2D_ts = types::list<regular_t, regular_t>;

        for (const grid_payload<material_slab_payload> &g_data :
             grids_data.grids) {
            // Start form finding the grid local frame and then build the grid
            if constexpr (DIM() == 2) {
                grid_reader_t::template deserialize<bounds2D_ts, binning2D_ts>(
                    g_data, det_builder);
            } else if constexpr (DIM() == 3) {
                using bounds3D_ts =
                    types::push_back<bounds2D_ts,
                                     closed<static_cast<label>(2u)>>;
                using binning3D_ts = types::push_back<binning2D_ts, regular_t>;

                grid_reader_t::template deserialize<bounds3D_ts, binning3D_ts>(
                    g_data, det_builder);
            } else {
                throw std::invalid_argument(
                    "No 1D material grid type defined in detray");
            }
        }
    }
};

}  // namespace detray
