/** Detray library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cstddef>

#include "detray/concepts/concept.hpp"

#ifdef DETRAY_HAVE_CONCEPTS
namespace detray::concepts {
template <typename T, std::size_t N>
concept transform = true;
}
#endif
