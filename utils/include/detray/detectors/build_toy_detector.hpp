/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/grid_builder.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_helper.hpp"
#include "detray/detectors/factories/barrel_generator.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/materials/mixture.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace detray {

/// Configure the toy detector
template <typename scalar_t>
struct toy_config {

    /// Default toy detector configuration
    toy_config()
        : m_n_brl_layers{4u},
          m_n_edc_layers{3u},
          m_use_material_maps{false},
          m_cyl_map_bins{20u, 20u},
          m_disc_map_bins{3u, 20u},
          m_mapped_material{
              mixture<scalar_t, silicon_tml<scalar_t, std::ratio<9, 10>>,
                      aluminium<scalar_t, std::ratio<1, 10>>>{}},
          m_thickness{1.5f * unit<scalar_t>::mm},
          m_cyl_mat_generator{detray::detail::generate_cyl_mat},
          m_disc_mat_generator{detray::detail::generate_disc_mat},
          m_do_check{true} {

        m_barrel_layer_radii = {
            19.f * unit<scalar_t>::mm, 32.f * unit<scalar_t>::mm,
            72.f * unit<scalar_t>::mm, 116.f * unit<scalar_t>::mm,
            172.f * unit<scalar_t>::mm};

        m_barrel_layer_sizes = {
            {0.f * unit<scalar_t>::mm, 27.f * unit<scalar_t>::mm},
            {27.f * unit<scalar_t>::mm, 38.f * unit<scalar_t>::mm},
            {64.f * unit<scalar_t>::mm, 80.f * unit<scalar_t>::mm},
            {108.f * unit<scalar_t>::mm, 124.f * unit<scalar_t>::mm},
            {164.f * unit<scalar_t>::mm, 180.f * unit<scalar_t>::mm}};

        m_barrel_binning = {
            {0u, 0u}, {16u, 14u}, {32u, 14u}, {52u, 14u}, {78u, 14u}};

        // Barrel module creation
        m_barrel_factory_cfg.half_length(500.f * unit<scalar_t>::mm)
            .module_bounds(
                {8.4f * unit<scalar_t>::mm, 36.f * unit<scalar_t>::mm})
            .tilt_phi(0.14f /*0.145*/)
            .radial_stagger(0.5f * unit<scalar_t>::mm /*2.f*/)
            .z_overlap(2.f * unit<scalar_t>::mm /*5.f*/);
    }

    /// No. of barrel layers the detector should be built with
    unsigned int m_n_brl_layers;
    /// No. of endcap layers (on either side) the detector should be built with
    unsigned int m_n_edc_layers;
    /// Do material maps on portals
    bool m_use_material_maps;
    /// Number of bins for material maps
    std::array<std::size_t, 2> m_cyl_map_bins;
    std::array<std::size_t, 2> m_disc_map_bins;
    /// Material to be filled into the maps
    material<scalar_t> m_mapped_material;
    /// Minimal thickness of the material slabs
    scalar_t m_thickness;
    /// Generate material along z bins for a cylinder material grid
    std::function<std::vector<material_slab<scalar_t>>(
        const std::array<scalar_t, 2u> &, const std::size_t, material<scalar_t>,
        const scalar_t)>
        m_cyl_mat_generator;
    /// Generate material along r bins for a disc material grid
    std::function<std::vector<material_slab<scalar_t>>(
        const std::array<scalar_t, 2u> &, const std::size_t, material<scalar_t>,
        const scalar_t)>
        m_disc_mat_generator;

    /// Radii at which to place the module layers (including beampipe)
    std::vector<scalar_t> m_barrel_layer_radii;
    /// Radii at which to place the layer portals
    std::vector<std::pair<scalar, scalar>> m_barrel_layer_sizes;
    /// Number of modules in phi and z for the barrel
    std::vector<std::pair<unsigned int, unsigned int>> m_barrel_binning;
    /// Config for the module generation
    barrel_generator_config<scalar_t> m_barrel_factory_cfg;
    /// Run detector consistency check after reading
    bool m_do_check;

    /// Setters
    /// @{
    constexpr toy_config &n_brl_layers(const unsigned int n) {
        m_n_brl_layers = n;
        return *this;
    }
    constexpr toy_config &n_edc_layers(const unsigned int n) {
        m_n_edc_layers = n;
        return *this;
    }
    constexpr toy_config &use_material_maps(const bool b) {
        m_use_material_maps = b;
        return *this;
    }
    constexpr toy_config &cyl_map_bins(const std::size_t n_rphi,
                                       const std::size_t n_z) {
        m_cyl_map_bins = {n_rphi, n_z};
        return *this;
    }
    constexpr toy_config &disc_map_bins(const std::size_t n_r,
                                        const std::size_t n_phi) {
        m_disc_map_bins = {n_r, n_phi};
        return *this;
    }
    constexpr toy_config &thickness(const scalar_t t) {
        assert(t > 0.f);
        m_thickness = t;
        return *this;
    }
    constexpr toy_config &mapped_material(const material<scalar_t> &mat) {
        m_mapped_material = mat;
        return *this;
    }
    constexpr toy_config &do_check(const bool check) {
        m_do_check = check;
        return *this;
    }
    /// @}

    /// Getters
    /// @{
    constexpr unsigned int n_brl_layers() const { return m_n_brl_layers; }
    constexpr unsigned int n_edc_layers() const { return m_n_edc_layers; }
    constexpr bool use_material_maps() const { return m_use_material_maps; }
    constexpr const std::array<std::size_t, 2> &cyl_map_bins() const {
        return m_cyl_map_bins;
    }
    constexpr const std::array<std::size_t, 2> &disc_map_bins() const {
        return m_disc_map_bins;
    }
    constexpr scalar_t thickness() const { return m_thickness; }
    auto barrel_mat_generator() const { return m_cyl_mat_generator; }
    auto edc_mat_generator() const { return m_disc_mat_generator; }
    constexpr material<scalar_t> mapped_material() const {
        return m_mapped_material;
    }
    constexpr const auto &barrel_layer_radii() const {
        return m_barrel_layer_radii;
    }
    constexpr const auto &barrel_layer_sizes() const {
        return m_barrel_layer_sizes;
    }
    constexpr const auto &barrel_layer_binning() const {
        return m_barrel_binning;
    }
    constexpr barrel_generator_config<scalar_t> &barrel_config() {
        return m_barrel_factory_cfg;
    }
    constexpr bool do_check() const { return m_do_check; }
    /// @}
};

/// Helper method for creating the barrel section.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param gctx geometry context
/// @param cfg config for the toy detector
/// @param beampipe_idx index of the beampipe outermost volume
/// @param n_layers the number of layers that should be built
/// @param factory_cfg config struct for module creation
template <typename detector_builder_t, typename scalar_t>
inline void add_barrel_detector(
    detector_builder_t &det_builder,
    typename detector_builder_t::detector_type::geometry_context &gctx,
    toy_config<scalar_t> &cfg,
    typename detector_builder_t::detector_type::name_map &names,
    dindex beampipe_idx) {

    using detector_t = typename detector_builder_t::detector_type;

    // Generate volume sizes in r, including gap volumes
    const auto &layer_sizes = cfg.barrel_layer_sizes();
    std::vector<std::pair<scalar_t, scalar_t>> vol_sizes{
        {layer_sizes[1].first, layer_sizes[1].second}};

    for (unsigned int i = 2u; i < cfg.n_brl_layers() + 1u; ++i) {
        vol_sizes.emplace_back(layer_sizes[i].first,
                               layer_sizes[i - 1u].second);
        vol_sizes.emplace_back(layer_sizes[i].first, layer_sizes[i].second);
    }

    // Alternate barrel module layers and gap volumes
    bool is_gap = true;
    for (unsigned int i = 0u; i < 2u * cfg.n_brl_layers() - 1u; ++i) {

        // New volume
        auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
        const dindex vol_idx{v_builder->vol_index()};

        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {

            /*det_helper.create_cyl_volume(cfg, det, resource, ctx,
                                         vol_sizes[i].first,
                                         vol_sizes[i].second, -brl_half_z,
                                         brl_half_z, volume_links_vec[i]);*/

            names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);
        } else {
            // Configure the module factory for this layer
            unsigned int j = (i + 2u) / 2u;
            const auto &barrel_binning = cfg.barrel_layer_binning().at(j);

            auto &barrel_cfg = cfg.barrel_config();
            barrel_cfg.binning(barrel_binning.first, barrel_binning.second);
            barrel_cfg.radius(cfg.barrel_layer_radii().at(j));

            auto mod_factory =
                std::make_shared<barrel_generator<detector_t, rectangle2D>>(
                    barrel_cfg);
            /*det_helper.create_cyl_volume(cfg, det, resource, ctx,
                                         vol_sizes[i].first,
                                         vol_sizes[i].second, -brl_half_z,
                                         brl_half_z, volume_links_vec[i]);*/

            v_builder->add_surfaces(mod_factory);

            names[vol_idx + 1u] = "barrel_" + std::to_string(vol_idx);

            // dindex_range sf_range = add_cylinder_grid(
            //     ctx, resource, det.volumes().back(), det, m_factory);
        }
    }
}

/// Builds a detray geometry that contains the innermost tml layers. The number
/// of barrel and endcap layers can be chosen, but all barrel layers should be
/// present when an endcap detector is built to have the barrel region radius
/// match the endcap diameter.
///
/// @param resource vecmem memory resource to use for container allocations
/// @param cfg toy detector configuration
///
/// @returns a complete detector object
template <typename scalar_t>
inline auto build_toy_detector(vecmem::memory_resource &resource,
                               toy_config<scalar_t> &cfg = {}) {

    using builder_t = detector_builder<toy_metadata, volume_builder>;
    using detector_t = typename builder_t::detector_type;
    // using material_id = typename detector_t::materials::id;
    // using nav_link_t = typename detector_t::surface_type::navigation_link;

    // Check config
    /*if (cfg.n_edc_layers() > edc_positions.size()) {
        throw std::invalid_argument(
            "ERROR: Too many endcap layers requested (max " +
            std::to_string(edc_positions.size()) + ")!");
    }*/
    if (cfg.n_brl_layers() > cfg.barrel_layer_radii().size() - 1u) {
        throw std::invalid_argument(
            "ERROR: Too many barrel layers requested (max " +
            std::to_string(cfg.barrel_layer_radii().size() - 1u) + ")!");
    }
    // the radius of the endcaps and  the barrel section need to match
    /*if (cfg.n_edc_layers() > 0 and
        math::abs(brl_lay_sizes[cfg.n_brl_layers()].second -
                  edc_config.outer_r) >
            std::numeric_limits<scalar>::epsilon()) {
        throw std::invalid_argument(
            "ERROR: Barrel and endcap radii do not match!");
    }*/

    // Toy detector builder
    builder_t det_builder;

    // Detector and volume names
    typename detector_t::name_map name_map = {{0u, "toy_detector"}};
    // Geometry context object
    typename detector_t::geometry_context gctx{};
    // Leaving world
    // constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    // Build the barrel section
    if (cfg.n_brl_layers() > 0u) {

        const dindex beampipe_idx{0u};

        add_barrel_detector(det_builder, gctx, cfg, name_map, beampipe_idx);
    }

    // Build and return the detector
    auto det = det_builder.build(resource);

    if (cfg.do_check()) {
        detray::detail::check_consistency(det);
    }

    return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
