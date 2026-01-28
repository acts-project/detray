/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <type_traits>

namespace detray {

/// Base class actor implementation
struct base_actor {
    /// Tag whether this is a composite type
    struct is_comp_actor : public std::false_type {};

    /// Defines the actors state. Hidden by actor implementations.
    struct state {};
};

}  // namespace detray
