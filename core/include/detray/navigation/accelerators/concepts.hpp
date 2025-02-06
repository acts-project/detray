// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s)
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/utils/grid/detail/concepts.hpp"

// System include(s)
#include <concepts>

namespace detray::concepts {

template <class accelerator_t>
concept surface_grid = concepts::grid<accelerator_t>&& std::same_as<
    typename accelerator_t::value_type,
    surface_descriptor<typename accelerator_t::value_type::mask_link,
                       typename accelerator_t::value_type::material_link,
                       typename accelerator_t::value_type::transform_link,
                       typename accelerator_t::value_type::navigation_link>>;

}  // namespace detray::concepts
