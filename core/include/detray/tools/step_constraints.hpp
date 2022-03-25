/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray definitions
#include <climits>
#include <type_traits>

#include "detray/definitions/qualifiers.hpp"

namespace detray {

// Add constraints for steppers
namespace step {

/// the types of constraints
/// from accuracy - this can vary up and down given a good step estimator
/// from actor    - this would be a typical navigation step
/// from aborter  - this would be a target condition
/// from user     - this is user given for what reason ever
enum constraint : std::size_t {
    accuracy = 0,
    actor = 1,
    aborter = 2,
    user = 3,
    all = 4
};

/// Direction in which the integration is performed
enum direction : int {
    e_forward = 1,
    e_unknown = 0,
    e_backward = -1,
};

}  // namespace step

/// Struct that represents unconstrained stepping
struct unconstrained_step {

    /// Register a new @param step_size constraint
    template <step::constraint type>
    DETRAY_HOST_DEVICE constexpr void register_constraint(
        scalar /*step_size*/) {}

    /// @returns the current step size constraint
    template <step::constraint type>
    DETRAY_HOST_DEVICE constexpr scalar size() {
        return std::numeric_limits<scalar>::max();
    }

    /// Remove constraints
    template <step::constraint type = step::constraint::all>
    DETRAY_HOST_DEVICE constexpr void release() {}
};

/// Struct that holds the current step size constraint, regardless of type.
struct single_constrainted_step {

    /// Register a new @param step_size constraint
    template <step::constraint type = step::constraint::actor>
    DETRAY_HOST_DEVICE void register_constraint(scalar step_size) {
        if (_direction == step::direction::e_unknown) {
            _direction = step_size > 0 ? step::direction::e_forward
                                       : step::direction::e_backward;
        }
        _constraint = std::min(_constraint, _direction * step_size);
    }

    /// @returns the current step size constraint
    template <step::constraint type>
    DETRAY_HOST_DEVICE scalar size() {
        return _direction * _constraint;
    }

    /// Remove all constraints
    template <step::constraint type>
    DETRAY_HOST_DEVICE void release() {
        _constraint = _direction * std::numeric_limits<scalar>::max();
        _direction = step::direction::e_unknown;
    }

    /// Current navigation direction. Can only be changed after step size
    /// release by first step size registration.
    step::direction _direction = step::direction::e_unknown;

    /// Current strongest step size constraint
    scalar _constraint = std::numeric_limits<scalar>::max();
};

/// Struct that can be configured with a number of different step sizes by other
/// actors and will then resolve the strictest one.
template <template <typename, std::size_t> class array_t = darray>
struct constrainted_step {

    /// Register a new @param step_size constraint
    template <step::constraint type,
              std::enable_if_t<not(type == step::constraint::all), bool> = true>
    DETRAY_HOST_DEVICE void register_constraint(scalar step_size) {
        if (_direction == step::direction::e_unknown) {
            _direction = step_size > 0 ? step::direction::e_forward
                                       : step::direction::e_backward;
        }
        _constraints[type] =
            std::min(_constraints[type], _direction * step_size);
    }

    /// @returns the current step size constraint for a given type or overall
    template <step::constraint type = step::constraint::all>
    DETRAY_HOST_DEVICE scalar size() {
        if constexpr (type == step::constraint::all) {
            return _direction * min();
        } else {
            return _direction * _constraints[type];
        }
    }

    template <step::constraint type = step::constraint::all>
    DETRAY_HOST_DEVICE void release(
        step::direction dir = step::direction::e_unknown) {
        if constexpr (type == step::constraint::all) {
            _constraints = {std::numeric_limits<scalar>::max()};
        } else {
            _constraints[type] = min();
        }
        _direction = dir;
    }

    /// @returns the strongest constraint
    inline scalar min() {
        scalar min_constr = std::numeric_limits<scalar>::max();
        min_constr =
            std::min(min_constr, _constraints[step::constraint::accuracy]);
        min_constr =
            std::min(min_constr, _constraints[step::constraint::actor]);
        min_constr =
            std::min(min_constr, _constraints[step::constraint::aborter]);
        return std::min(min_constr, _constraints[step::constraint::user]);
    }

    /// Current navigation direction. Can only be changed after step size
    /// release by first step size registration.
    step::direction _direction = step::direction::e_unknown;

    /// Current strongets step size constraint
    array_t<scalar, 4> _constraints = {std::numeric_limits<scalar>::max()};
};

}  // namespace detray