/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/propagator/propagation_config.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray::test {

/// Test configuration type
template <concepts::scalar scalar_t>
struct configuration {
    /// General testing
    /// @{
    /// Tolerance to compare two floating point values
    scalar_t m_tolerance{std::numeric_limits<scalar_t>::epsilon()};
    /// Shorthand for infinity
    scalar_t inf{std::numeric_limits<scalar_t>::infinity()};
    /// Shorthand for the floating point epsilon
    scalar_t epsilon{std::numeric_limits<scalar_t>::epsilon()};
    /// @}

    /// Propagation
    /// @{
    propagation::config m_prop_cfg{};
    /// @}

    /// Setters
    /// @{
    configuration& tol(scalar_t t) {
        m_tolerance = t;
        return *this;
    }
    /// @}

    /// Getters
    /// @{
    scalar_t tol() const { return m_tolerance; }
    propagation::config& propagation() { return m_prop_cfg; }
    const propagation::config& propagation() const { return m_prop_cfg; }
    /// @}

    /// Print configuration
    std::ostream& operator<<(std::ostream& os) {
        os << " -> test tolerance:  \t " << tol() << std::endl;
        os << " -> trk path limit:  \t " << propagation().stepping.path_limit
           << std::endl;
        os << " -> overstepping tol:\t "
           << propagation().navigation.overstep_tolerance << std::endl;
        os << " -> step constraint:  \t "
           << propagation().stepping.step_constraint << std::endl;
        os << std::endl;

        return os;
    }
};

}  // namespace detray::test
