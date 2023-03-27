/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray::test {

/// Base type for linear algebra benchmarks with google benchmark
struct fixture_base {
    /// Useful typedefs
    /// @{
    using scalar_type = detray::scalar;
    using point2_type = __plugin::point2<scalar>;
    using point3_type = __plugin::point3<scalar>;
    using vector3_type = __plugin::vector3<scalar>;
    using transform3_type = __plugin::transform3<detray::scalar>;
    using matrix_operator = typename transform3_type::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    /// @}

    /// Local configuration type
    struct configuration {
        /// General testing
        /// @{
        /// Tolerance to compare two floating point values
        scalar_type m_tolerance{std::numeric_limits<scalar_type>::epsilon()};
        /// Shorthand for infinity
        scalar_type inf{std::numeric_limits<scalar_type>::infinity()};
        /// Shorthand for the floating point epsilon
        scalar_type epsilon{std::numeric_limits<scalar_type>::epsilon()};
        /// @}

        /// Propagation
        /// @{
        /// Maximum total pathlength for a track
        scalar_type m_path_limit{5.f * unit<scalar>::cm};
        scalar_type m_overstep_tolerance{-5.f * unit<scalar>::um};
        scalar_type m_step_constraint{std::numeric_limits<scalar_type>::max()};
        /// @}

        /// Setters
        /// @{
        configuration& tol(scalar t) {
            m_tolerance = t;
            return *this;
        }
        configuration& path_limit(scalar lim) {
            assert(lim > 0.f);
            m_tolerance = lim;
            return *this;
        }
        configuration& overstepping_tolerance(scalar t) {
            assert(t <= 0.f);
            m_tolerance = t;
            return *this;
        }
        configuration& step_containt(scalar constr) {
            assert(constr > 0.f);
            m_tolerance = constr;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        scalar_type tol() const { return m_tolerance; }
        scalar_type path_limit() const { return m_path_limit; }
        scalar_type overstepping_tolerance() const {
            return m_overstep_tolerance;
        }
        scalar_type step_containt() const { return m_step_constraint; }
        /// @}

        /// Print configuration
        friend std::ostream& operator<<(std::ostream& os,
                                        const configuration& c);
    };

    /// @returns the benchmark name
    std::string name() const { return "detray_test"; };
};

std::ostream& operator<<(std::ostream& os,
                         const fixture_base::configuration& cfg) {
    os << " -> test tolerance:  \t " << cfg.tol() << std::endl;
    os << " -> trk path limit:  \t " << cfg.path_limit() << std::endl;
    os << " -> overstepping tol:\t " << cfg.overstepping_tolerance()
       << std::endl;
    os << " -> step constraint:  \t " << cfg.step_containt() << std::endl;
    os << std::endl;
    return os;
}

}  // namespace detray::test