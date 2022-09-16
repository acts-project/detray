/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "detray/utils/ranges/enumerate.hpp"
#include "detray/utils/ranges/iota.hpp"

namespace detray::views {

/// Shorthand to test a type for being iterable.
template <class T>
using iota = detray::ranges::iota_view<T>;

}  // namespace detray::views
