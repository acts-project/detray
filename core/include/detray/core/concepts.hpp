// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s)
#include "detray/core/detail/type_traits.hpp"

namespace detray::concepts {

/// Check for the the presence of any type of grids in a detector definition
template <class detector_t>
concept has_grids = detail::contains_grids_v<typename detector_t::accel> ||
                    detail::contains_grids_v<typename detector_t::materials>;

/// Check for the the presence of surface grids in a detector definition
template <class detector_t>
concept has_surface_grids =
    detail::contains_surface_grids_v<typename detector_t::accel>;

/// Check for the the presence of material slabs in a detector definition
template <class detector_t>
concept has_material_slabs =
    detail::contains_material_slabs_v<typename detector_t::materials>;

/// Check for the the presence of material rods in a detector definition
template <class detector_t>
concept has_material_rods =
    detail::contains_material_rods_v<typename detector_t::materials>;

/// Check for the the presence of homogeneous material types in a detector
/// definition
template <class detector_t>
concept has_homogeneous_material =
    detail::contains_homogeneous_material_v<typename detector_t::materials>;

/// Check for the the presence of material maps in a detector definition
template <class detector_t>
concept has_material_maps =
    detail::contains_material_maps_v<typename detector_t::materials>;

}  // namespace detray::concepts
