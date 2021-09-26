/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/populator.hpp"
#include "grids/serializer2.hpp"
#include "utils/indexing.hpp"

#pragma once

namespace detray {

static constexpr int n_points = 3;

using host_replace = host_replace_populator<test::point3>;
using device_replace = device_replace_populator<test::point3>;

using host_complete = host_complete_populator<n_points, false, test::point3>;
using device_complete =
    device_complete_populator<n_points, false, test::point3>;

using host_attach = host_attach_populator<false, test::point3>;
using device_attach = device_attach_populator<false, test::point3>;

}  // namespace detray
