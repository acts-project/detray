/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/propagator.hpp"

// Test include(s)
#include "detray/test/common/types.hpp"

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
    /// Linear algebra typedefs
    /// @{
    using algebra_type = ALGEBRA_PLUGIN<test::scalar>;
    using scalar_type = dscalar<algebra_type>;
    using point2_type = dpoint2D<algebra_type>;
    using point3_type = dpoint3D<algebra_type>;
    using vector3_type = dvector3D<algebra_type>;
    using transform3_type = dtransform3D<algebra_type>;
    /// @}

    /// Local configuration type
    struct configuration {
        /// General testing
        /// @{
        /// Tolerance to compare two floating point values
        float m_tolerance{std::numeric_limits<float>::epsilon()};
        /// Shorthand for infinity
        float inf{std::numeric_limits<float>::infinity()};
        /// Shorthand for the floating point epsilon
        float epsilon{std::numeric_limits<float>::epsilon()};
        /// @}

        /// Propagation
        /// @{
        propagation::config m_prop_cfg{};
        /// @}

        /// Setters
        /// @{
        configuration& tol(float t) {
            m_tolerance = t;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        float tol() const { return m_tolerance; }
        propagation::config& propagation() { return m_prop_cfg; }
        const propagation::config& propagation() const { return m_prop_cfg; }
        /// @}

        /// Print configuration
        std::ostream& operator<<(std::ostream& os) {
            os << " -> test tolerance:  \t " << tol() << std::endl;
            os << " -> trk path limit:  \t "
               << propagation().stepping.path_limit << std::endl;
            os << " -> overstepping tol:\t "
               << propagation().navigation.overstep_tolerance << std::endl;
            os << " -> step constraint:  \t "
               << propagation().stepping.step_constraint << std::endl;
            os << std::endl;

            return os;
        }
    };

    /// Constructor
    fixture_base(const configuration& cfg = {})
        : tolerance{cfg.tol()},
          inf{cfg.inf},
          epsilon{cfg.epsilon},
          path_limit{cfg.propagation().stepping.path_limit},
          overstep_tolerance{cfg.propagation().navigation.overstep_tolerance},
          step_constraint{cfg.propagation().stepping.step_constraint} {}

    /// @returns the benchmark name
    std::string name() const { return "detray_test"; };

    protected:
    float tolerance{}, inf{}, epsilon{}, path_limit{}, overstep_tolerance{},
        step_constraint{};

    //static void SetUpTestSuite() {}
    //static void TearDownTestSuite() {}
    //virtual void SetUp() override {}
    //virtual void TearDown() override {}
};

}  // namespace detray::test
