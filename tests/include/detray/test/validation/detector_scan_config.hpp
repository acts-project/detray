/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray plugin include(s)
#include "detray/plugins/svgtools/styling/styling.hpp"

// Detray test include(s)
#include "detray/test/common/test_configuration.hpp"

// System include(s)
#include <limits>
#include <memory>
#include <string>

namespace detray::test {

/// @brief Configuration for a detector scan test.
template <typename track_generator_t, concepts::algebra algebra_t>
struct detector_scan_config
    : public detray::test::configuration<dscalar<algebra_t>> {

    using scalar_type = dscalar<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using base_type = detray::test::configuration<scalar_type>;
    using trk_gen_config_t = typename track_generator_t::configuration;

    /// Name of the test
    std::string m_name{"detector_scan"};
    /// Name of the input file, containing the complete ray scan traces
    std::string m_intersection_file{"truth_intersections"};
    std::string m_track_param_file{"truth_trk_parameters"};
    /// Mask tolerance for the intersectors
    darray<scalar_type, 2> m_mask_tol{
        std::numeric_limits<scalar_type>::epsilon(),
        std::numeric_limits<scalar_type>::epsilon()};
    /// B-field vector for helix
    vector3_type m_B{0.f * unit<scalar_type>::T, 0.f * unit<scalar_type>::T,
                     2.f * unit<scalar_type>::T};
    /// Configuration of the ray generator
    trk_gen_config_t m_trk_gen_cfg{};
    /// Write intersection points for plotting
    bool m_write_inters{false};
    /// Visualization style to be applied to the svgs
    detray::svgtools::styling::style m_style =
        detray::svgtools::styling::tableau_colorblind::style;

    /// Getters
    /// @{
    const std::string &name() const { return m_name; }
    const std::string &intersection_file() const { return m_intersection_file; }
    const std::string &track_param_file() const { return m_track_param_file; }
    darray<scalar_type, 2> mask_tolerance() const { return m_mask_tol; }
    const vector3_type &B_vector() { return m_B; }
    bool write_intersections() const { return m_write_inters; }
    trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
    const trk_gen_config_t &track_generator() const { return m_trk_gen_cfg; }
    const auto &svg_style() const { return m_style; }
    /// @}

    /// Setters
    /// @{
    detector_scan_config &name(const std::string_view n) {
        m_name = n;
        return *this;
    }
    detector_scan_config &intersection_file(const std::string_view f) {
        m_intersection_file = f;
        return *this;
    }
    detector_scan_config &track_param_file(const std::string_view f) {
        m_track_param_file = f;
        return *this;
    }
    detector_scan_config &mask_tolerance(const darray<scalar_type, 2> tol) {
        m_mask_tol = tol;
        return *this;
    }
    detector_scan_config &B_vector(const vector3_type &B) {
        m_B = B;
        return *this;
    }
    detector_scan_config &B_vector(const scalar_type B0, const scalar_type B1,
                                   const scalar_type B2) {
        m_B = vector3_type{B0, B1, B2};
        return *this;
    }
    detector_scan_config &write_intersections(const bool do_write) {
        m_write_inters = do_write;
        return *this;
    }
    /// @}
};

}  // namespace detray::test
