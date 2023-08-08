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
        /// Maximum total pathlength for a track
        scalar m_path_limit{5.f * unit<scalar>::cm};
        scalar m_overstep_tolerance{-5.f * unit<scalar>::um};
        scalar m_step_constraint{std::numeric_limits<scalar>::max()};
        /// @}

        /// Setters
        /// @{
        configuration& tol(scalar t) {
            m_tolerance = t;
            return *this;
        }
        configuration& path_limit(scalar lim) {
            assert(lim > 0.f);
            m_path_limit = lim;
            return *this;
        }
        configuration& overstepping_tolerance(scalar t) {
            assert(t <= 0.f);
            m_overstep_tolerance = t;
            return *this;
        }
        configuration& step_constraint(scalar constr) {
            assert(constr > 0.f);
            m_step_constraint = constr;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        scalar tol() const { return m_tolerance; }
        scalar path_limit() const { return m_path_limit; }
        scalar overstepping_tolerance() const { return m_overstep_tolerance; }
        scalar step_constraint() const { return m_step_constraint; }
        /// @}

        /// Print configuration
        std::ostream& operator<<(std::ostream& os) {
            os << " -> test tolerance:  \t " << tol() << std::endl;
            os << " -> trk path limit:  \t " << path_limit() << std::endl;
            os << " -> overstepping tol:\t " << overstepping_tolerance()
               << std::endl;
            os << " -> step constraint:  \t " << step_constraint() << std::endl;
            os << std::endl;

            return os;
        }
    };

    /// Constructor
    fixture_base(const configuration& cfg = {})
        : tolerance{cfg.tol()},
          inf{cfg.inf},
          epsilon{cfg.epsilon},
          path_limit{cfg.path_limit()},
          overstep_tolerance{cfg.overstepping_tolerance()},
          step_constraint{cfg.step_constraint()} {}

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
