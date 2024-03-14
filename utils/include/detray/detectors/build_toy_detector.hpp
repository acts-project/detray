/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/cylinder_portal_generator.hpp"
#include "detray/builders/detail/portal_accessor.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/grid_builder.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_helper.hpp"
#include "detray/detectors/factories/barrel_generator.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
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

namespace detail {

/// Add the cylinder and disc portals for a volume from explicit parameters
///
/// @param v_builder the volume builder to add the portals to
/// @param names name map for volumes of the detector under construction
/// @param inner_r inner volume radius
/// @param outer_r outer volume radius
/// @param h_z half length of the cylinder
/// @param link_north portal volume link of the outer cylinder
/// @param link_south portal volume link of the inner cylinder
/// @param link_east portal volume link of the left disc
/// @param link_west portal volume link of the right disc
template <typename detector_t, typename scalar_t>
void add_gap_portals(volume_builder_interface<detector_t> *v_builder,
                     typename detector_t::name_map &names,
                     const scalar_t inner_r, const scalar_t outer_r,
                     const scalar_t h_z, const dindex link_north,
                     const dindex link_south, const dindex link_east,
                     const dindex link_west) {

    using transform3_t = typename detector_t::transform3;
    using point3_t = typename detector_t::point3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    const transform3_t identity{};
    const dindex vol_idx{v_builder->vol_index()};

    using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
    using disc_factory_t = surface_factory<detector_t, ring2D>;

    auto pt_cyl_factory = std::make_shared<cyl_factory_t>();
    auto pt_disc_factory = std::make_shared<disc_factory_t>();

    // Inner cylinder portal
    pt_cyl_factory->push_back({surface_id::e_portal, identity,
                               static_cast<nav_link_t>(link_south),
                               std::vector<scalar_t>{inner_r, -h_z, h_z}});
    // Outer cylinder portal
    pt_cyl_factory->push_back({surface_id::e_portal, identity,
                               static_cast<nav_link_t>(link_north),
                               std::vector<scalar_t>{outer_r, -h_z, h_z}});

    // Left disc portal
    pt_disc_factory->push_back({surface_id::e_portal,
                                transform3_t{point3_t{0.f, 0.f, -h_z}},
                                static_cast<nav_link_t>(link_east),
                                std::vector<scalar_t>{inner_r, outer_r}});
    // Right disc portal
    pt_disc_factory->push_back({surface_id::e_portal,
                                transform3_t{point3_t{0.f, 0.f, h_z}},
                                static_cast<nav_link_t>(link_west),
                                std::vector<scalar_t>{inner_r, outer_r}});

    v_builder->add_surfaces(pt_cyl_factory);
    v_builder->add_surfaces(pt_disc_factory);

    names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);
}

/// Helper method for creating the barrel section.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param gctx geometry context
/// @param cfg config for the toy detector
/// @param names name map for volumes of the detector under construction
/// @param beampipe_idx index of the beampipe outermost volume
template <typename detector_builder_t, typename scalar_t>
inline void add_barrel_detector(
    detector_builder_t &det_builder,
    typename detector_builder_t::detector_type::geometry_context &gctx,
    toy_config<scalar_t> &cfg,
    typename detector_builder_t::detector_type::name_map &names,
    dindex beampipe_idx) {

    using detector_t = typename detector_builder_t::detector_type;
    using transform3_t = typename detector_t::transform3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
    const transform3_t identity{};

    const scalar h_z{cfg.barrel_config().half_length()};
    // Set the inner radius of the first gat to the radius of the beampipe vol.
    scalar_t gap_inner_r{cfg.barrel_layer_radii().at(0)};

    typename cylinder_portal_generator<detector_t>::boundaries vol_bounds{};

    // Alternate barrel module layers and gap volumes
    bool is_gap = true;
    for (unsigned int i = 0u; i < 2u * cfg.n_brl_layers(); ++i) {

        // New volume
        auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
        const dindex vol_idx{v_builder->vol_index()};
        // The barrel volumes are centered at the origin
        v_builder->add_volume_placement(identity);

        auto link_north{vol_idx + 1u};
        auto link_south{vol_idx - 1u};
        auto link_east{end_of_world};
        auto link_west{end_of_world};

        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {

            // The first time a gap is built, it needs to link to the beampipe
            link_south = (i == 2u) ? beampipe_idx : link_south;

            detail::add_gap_portals(v_builder, names, gap_inner_r,
                                    vol_bounds.inner_radius, h_z, link_north,
                                    link_south, link_east, link_west);

            // Set the inner radius for the next gap volume
            gap_inner_r = vol_bounds.outer_radius;

        } else {
            // Configure the module factory for this layer
            auto &barrel_cfg = cfg.barrel_config();

            const unsigned int j{(i + 2u) / 2u};
            barrel_cfg.binning(cfg.barrel_layer_binning().at(j));
            barrel_cfg.radius(cfg.barrel_layer_radii().at(j));

            // Add a layer of module surfaces
            auto module_factory =
                std::make_shared<barrel_generator<detector_t, rectangle2D>>(
                    barrel_cfg);

            // Configure the portal factory
            cylinder_portal_config<scalar_t> portal_cfg{};
            portal_cfg.autofit(true).fixed_half_length(h_z);
            // Link the volume portals to its neighbors
            portal_cfg.link_north(link_north).link_south(link_south);
            portal_cfg.link_east(link_east).link_west(link_west);

            // Add cylinder and disc portals
            auto portal_factory =
                std::make_shared<cylinder_portal_generator<detector_t>>(
                    portal_cfg);

            v_builder->add_surfaces(module_factory);
            v_builder->add_surfaces(portal_factory);

            // Set the new current boundaries, to cunstruct the next gap
            vol_bounds = portal_factory->volume_boundaries();
            std::cout << "layer: " << vol_bounds.inner_radius << ", "
                      << vol_bounds.outer_radius << std::endl;

            names[vol_idx + 1u] = "barrel_" + std::to_string(vol_idx);

            // dindex_range sf_range = add_cylinder_grid(
            //     ctx, resource, det.volumes().back(), det, m_factory);
        }
    }

    // Add a final gap volume to get to the full barrel radius
    auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
    const dindex vol_idx{v_builder->vol_index()};
    v_builder->add_volume_placement(identity);

    detail::add_gap_portals(v_builder, names, vol_bounds.outer_radius,
                            180.f * unit<scalar_t>::mm, h_z, vol_idx - 1u,
                            end_of_world, end_of_world, end_of_world);
}

}  // namespace detail

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

        detail::add_barrel_detector(det_builder, gctx, cfg, name_map,
                                    beampipe_idx);
    }

    // Build and return the detector
    auto det = det_builder.build(resource);

    if (cfg.do_check()) {
        detray::detail::check_consistency(det);
    }

    return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
