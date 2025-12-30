/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace algebra::boolean {

/// Utilities for single booleans: default case
/// @{
constexpr bool any_of(bool b) {
  return b;
}
constexpr bool all_of(bool b) {
  return b;
}
constexpr bool none_of(bool b) {
  return !b;
}
/// @}

}  // namespace algebra::boolean
