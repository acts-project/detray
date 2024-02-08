/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/grid_reader.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/grid_builder.hpp"

// System include(s)
#include <string>

namespace detray {

/// @brief Abstract base class for surface grid readers
template <class detector_t,
          typename CAP = std::integral_constant<std::size_t, 9>,
          typename DIM = std::integral_constant<std::size_t, 2>>
class surface_grid_reader
    : public detail::grid_reader<detector_t, typename detector_t::surface_type,
                                 grid_builder, CAP, DIM> {

    using grid_reader_t =
        detail::grid_reader<detector_t, typename detector_t::surface_type,
                            grid_builder, CAP, DIM>;
    using base_type = grid_reader_t;

    protected:
    /// Tag the reader as "surface_grids"
    inline static const std::string tag = "surface_grids";

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
        const detector_grids_payload<std::size_t, io::detail::acc_type>
            &grids_data) {

        grid_reader_t::deserialize(det_builder, grids_data);
    }
};

}  // namespace detray
