/** Detray library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/concepts/concept.hpp"

#ifdef DETRAY_HAVE_CONCEPTS
namespace detray::concepts {
template <typename T>
concept axis = true;
}
#endif
