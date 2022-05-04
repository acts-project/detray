/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/materials/material.hpp"

namespace detray {

/// Homegenous material
template <int ID>
struct material<e_homogeneous, ID> {};

}  // namespace detray