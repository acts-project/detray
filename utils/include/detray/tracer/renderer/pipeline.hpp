/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/propagator/actor_chain.hpp"
#include "detray/tracer/renderer/colorizer.hpp"
#include "detray/tracer/renderer/intersector.hpp"
#include "detray/tracer/renderer/single_shape.hpp"

namespace detray {

template <template <typename...> class tuple_t = dtuple, typename... shaders_t>
using rendering_pipeline = actor_chain<tuple_t, shaders_t...>;

}  // namespace detray
