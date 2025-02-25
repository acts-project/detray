/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/navigation/direct_navigator.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/utils/detectors/build_toy_detector.hpp"
#include "detray/test/utils/detectors/build_wire_chamber.hpp"
#include "detray/test/utils/inspectors.hpp"

// Test include(s)
#include "detray/test/utils/types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GoogleTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <map>

namespace detray {

constexpr std::size_t cache_size{navigation::default_cache_size};

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {
    using context_t = typename navigation_t::detector_type::geometry_context;
    stepping_t _stepping;
    navigation_t _navigation;
    context_t _context{};
};

}  // namespace detray