/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <type_traits>

namespace detray::n_axis {

/// axis bounds names.
///
/// open: adds additional over-/undeflow bins.
/// closed: over-/underflow values are mapped into the bin range.
/// circular: over-/underflow values wrap around.
enum class bounds {
    e_open = 0,
    e_closed = 1,
    e_circular = 2,
};

/// axis coordinate names. Used to get a specific axes from an axes collection.
///
/// x, y, z: cartesian coordinate axes.
/// r, phi: polar coordinate axes.
/// rphi, cyl_z: 2D cylinder axes (3D cylinder uses r, phi, z).
enum class label {
    e_x = 0,
    e_y = 1,
    e_z = 2,
    e_r = 0,
    e_phi = 1,
    e_rphi = 0,
    e_cyl_z = 1,
};

/// axis binning type names.
///
/// regular: same sized bins along the axis.
/// irregular: every bin can have a different size along the axis.
enum class binning {
    e_regular = 0,
    e_irregular = 1,
};

/// Determine axis bounds as either 'open' or 'closed' for non-circular axes.
/// Used in the shape structs.
/// @{
template <n_axis::label L>
struct open;

template <n_axis::label L>
struct closed;

template <n_axis::bounds s, n_axis::label axis_label>
using bounds_t =
    std::conditional_t<s == n_axis::bounds::e_open, n_axis::open<axis_label>,
                       n_axis::closed<axis_label>>;
/// @}

}  // namespace detray::n_axis