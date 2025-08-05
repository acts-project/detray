/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/propagator.hpp"

// Detray test include(s)
#include "detray/test/common/test_configuration.hpp"
#include "detray/test/framework/types.hpp"

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
    using algebra_t = test::algebra;
    using scalar = test::scalar;
    using point2 = test::point2;
    using point3 = test::point3;
    using vector3 = test::vector3;
    using transform3 = test::transform3;
    /// @}

    using configuration = detray::test::configuration<scalar>;

    /// Constructor
    explicit fixture_base(const configuration& cfg = {})
        : tolerance{cfg.tol()},
          inf{cfg.inf},
          epsilon{cfg.epsilon},
          path_limit{cfg.propagation().stepping.path_limit},
          overstep_tolerance{cfg.propagation().navigation.overstep_tolerance},
          step_constraint{cfg.propagation().stepping.step_constraint} {}

    /// @returns the benchmark name
    std::string name() const { return "detray_test"; };

    protected:
    static void SetUpTestSuite() { /* Do nothing */ }
    static void TearDownTestSuite() { /* Do nothing */ }
    void SetUp() override { /* Do nothing */ }
    void TearDown() override { /* Do nothing */ }

    private:
    scalar tolerance{};
    scalar inf{};
    scalar epsilon{};
    scalar path_limit{};
    scalar overstep_tolerance{};
    scalar step_constraint{};
};

}  // namespace detray::test
