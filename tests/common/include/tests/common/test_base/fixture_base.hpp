/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/propagator.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray::test {

/// Base type for test fixtures with google test
template <typename scope = ::testing::Test>
class fixture_base : public scope {
    public:
    /// Useful typedefs
    /// @{
    using scalar = detray::scalar;
    using point2 = __plugin::point2<scalar>;
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using transform3 = __plugin::transform3<detray::scalar>;
    /// @}

    /// Local configuration type
    struct configuration {
        /// General testing
        /// @{
        /// Tolerance to compare two floating point values
        scalar m_tolerance{std::numeric_limits<scalar>::epsilon()};
        /// Shorthand for infinity
        scalar inf{std::numeric_limits<scalar>::infinity()};
        /// Shorthand for the floating point epsilon
        scalar epsilon{std::numeric_limits<scalar>::epsilon()};
        /// @}

        /// Propagation
        /// @{
        propagation::config<scalar> m_prop_cfg{};
        /// @}

        /// Setters
        /// @{
        configuration& tol(scalar t) {
            m_tolerance = t;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        scalar tol() const { return m_tolerance; }
        propagation::config<scalar>& propagation() { return m_prop_cfg; }
        const propagation::config<scalar>& propagation() const {
            return m_prop_cfg;
        }
        /// @}

        /// Print configuration
        std::ostream& operator<<(std::ostream& os) {
            os << " -> test tolerance:  \t " << tol() << std::endl;
            os << " -> trk path limit:  \t " << propagation().path_limit
               << std::endl;
            os << " -> overstepping tol:\t " << propagation().overstep_tolerance
               << std::endl;
            os << " -> step constraint:  \t " << propagation().step_constraint
               << std::endl;
            os << std::endl;

            return os;
        }
    };

    /// Constructor
    fixture_base(const configuration& cfg = {})
        : tolerance{cfg.tol()},
          inf{cfg.inf},
          epsilon{cfg.epsilon},
          path_limit{cfg.propagation().path_limit},
          overstep_tolerance{cfg.propagation().overstep_tolerance},
          step_constraint{cfg.propagation().step_constraint} {}

    /// @returns the benchmark name
    std::string name() const { return "detray_test"; };

    protected:
    scalar tolerance{}, inf{}, epsilon{}, path_limit{}, overstep_tolerance{},
        step_constraint{};

    static void SetUpTestSuite() {}
    static void TearDownTestSuite() {}
    virtual void SetUp() override {}
    virtual void TearDown() override {}
};

}  // namespace detray::test
