/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/units.hpp"

// System include(s).
#include <limits>

namespace detray::stepping {

enum class id {
    // False for non-charged tracks
    e_linear = 0,
    // True for charged tracks
    e_rk = 1,
};

template <typename scalar_t>
struct config {
    /// Minimum step size
    scalar_t min_stepsize{1e-4f * unit<scalar_t>::mm};
    /// Runge-Kutta numeric error tolerance
    scalar_t rk_error_tol{1e-4f * unit<scalar_t>::mm};
    /// Step size constraint
    scalar_t step_constraint{std::numeric_limits<scalar_t>::max()};
    /// Maximal path length of track
    scalar_t path_limit{5.f * unit<scalar_t>::m};
    /// Maximum number of Runge-Kutta step trials
    std::size_t max_rk_updates{10000u};
    /// Use mean energy loss (Bethe)
    /// if false, most probable energy loss (Landau) will be used
    bool use_mean_loss{true};
    /// Use eloss gradient in error propagation
    bool use_eloss_gradient{false};
    /// Use b field gradient in error propagation
    bool use_field_gradient{false};
};

}  // namespace detray::stepping
