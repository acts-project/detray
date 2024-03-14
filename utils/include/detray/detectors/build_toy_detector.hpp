/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/cylinder_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/grid_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_generator.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_helper.hpp"
#include "detray/detectors/factories/barrel_generator.hpp"
#include "detray/detectors/factories/endcap_generator.hpp"
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
#include <type_traits>
#include <utility>

namespace detray {

/// Configure the toy detector
template <typename scalar_t>
struct toy_config {

    /// Default toy detector configuration
    toy_config() {
        // Barrel module creator
        m_barrel_factory_cfg.half_length(500.f * unit<scalar_t>::mm)
            .module_bounds(
                {8.4f * unit<scalar_t>::mm, 36.f * unit<scalar_t>::mm})
            .tilt_phi(0.14f /*0.145*/)
            .radial_stagger(0.5f * unit<scalar_t>::mm /*2.f*/)
            .z_overlap(2.f * unit<scalar_t>::mm /*5.f*/);

        // Endcap module creator
        m_endcap_factory_cfg.inner_radius(m_beampipe_volume_radius)
            .outer_radius(m_outer_radius)
            .module_bounds(
                {{3.f * unit<scalar_t>::mm, 9.5f * unit<scalar_t>::mm,
                  39.f * unit<scalar_t>::mm},
                 {6.f * unit<scalar_t>::mm, 10.f * unit<scalar_t>::mm,
                  39.f * unit<scalar_t>::mm}})
            .ring_stagger(2.f * unit<scalar_t>::mm)
            .phi_stagger({4.f * unit<scalar_t>::mm, 4.f * unit<scalar_t>::mm})
            .phi_sub_stagger(
                {0.5f * unit<scalar_t>::mm, 0.5f * unit<scalar_t>::mm})
            .module_tilt({0.f, 0.f})
            .binning({40u, 68u});

        // Configure the material generation
        m_material_config.sensitive_material(silicon_tml<scalar_t>())
            .passive_material(beryllium_tml<scalar_t>())  // < beampipe
            .portal_material(vacuum<scalar_t>())
            .thickness(1.5f * unit<scalar_t>::mm);
    }

    /// No. of barrel layers the detector should be built with
    unsigned int m_n_brl_layers{4u};
    /// No. of endcap layers (on either side) the detector should be built with
    unsigned int m_n_edc_layers{3u};
    /// Total outer radius of the pixel subdetector
    scalar_t m_outer_radius{180.f * unit<scalar_t>::mm};
    // Radius of the innermost volume that contains the beampipe
    scalar_t m_beampipe_volume_radius{25.f * unit<scalar_t>::mm};
    // Envelope around the modules used by the cylinder portal generator
    scalar_t m_portal_envelope{0.5f * unit<scalar_t>::mm};
    /// Configuration for the homogeneous material generator
    hom_material_config<scalar_t> m_material_config{};
    /// Put material maps on portals or use homogenous material on modules
    bool m_use_material_maps{false};
    /// Number of bins for material maps
    std::array<std::size_t, 2> m_cyl_map_bins{20u, 20u};
    std::array<std::size_t, 2> m_disc_map_bins{3u, 20u};
    /// Material to be filled into the maps
    material<scalar_t> m_mapped_material{
        mixture<scalar_t, silicon_tml<scalar_t, std::ratio<9, 10>>,
                aluminium<scalar_t, std::ratio<1, 10>>>{}};
    /// Minimal thickness of the material slabs in the material maps
    scalar_t m_thickness{0.15f * unit<scalar_t>::mm};
    /// Thickness of the beampipe material
    scalar_t m_beampipe_mat_thickness{0.8f * unit<scalar_t>::mm};
    /// Thickness of the material slabs in the homogeneous material description
    scalar_t m_module_mat_thickness{1.5f * unit<scalar_t>::mm};
    /// Generate material along z bins for a cylinder material grid
    std::function<std::vector<material_slab<scalar_t>>(
        const std::array<scalar_t, 2u> &, const std::size_t, material<scalar_t>,
        const scalar_t)>
        m_cyl_mat_generator{detray::detail::generate_cyl_mat};
    /// Generate material along r bins for a disc material grid
    std::function<std::vector<material_slab<scalar_t>>(
        const std::array<scalar_t, 2u> &, const std::size_t, material<scalar_t>,
        const scalar_t)>
        m_disc_mat_generator{detray::detail::generate_disc_mat};
    /// Radii at which to place the barrel module layers (including beampipe)
    std::vector<scalar_t> m_barrel_layer_radii = {
        19.f * unit<scalar_t>::mm, 32.f * unit<scalar_t>::mm,
        72.f * unit<scalar_t>::mm, 116.f * unit<scalar_t>::mm,
        172.f * unit<scalar_t>::mm};
    /// Number of modules in phi and z for the barrel
    std::vector<std::pair<unsigned int, unsigned int>> m_barrel_binning = {
        {0u, 0u}, {16u, 14u}, {32u, 14u}, {52u, 14u}, {78u, 14u}};
    /// Positions at which to place the endcap module layers on either side
    std::vector<scalar_t> m_endcap_layer_positions = {
        600.f * unit<scalar_t>::mm,  700.f * unit<scalar_t>::mm,
        820.f * unit<scalar_t>::mm,  960.f * unit<scalar_t>::mm,
        1100.f * unit<scalar_t>::mm, 1300.f * unit<scalar_t>::mm,
        1500.f * unit<scalar_t>::mm};
    /// Config for the module generation (barrel)
    barrel_generator_config<scalar_t> m_barrel_factory_cfg{};
    /// Config for the module generation (endcaps)
    endcap_generator_config<scalar_t> m_endcap_factory_cfg{};
    /// Run detector consistency check after reading
    bool m_do_check{true};

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
    constexpr toy_config &envelope(const scalar_t env) {
        m_portal_envelope = env;
        return *this;
    }
    constexpr toy_config &use_material_maps(const bool b) {
        m_use_material_maps = b;
        return *this;
    }
    constexpr toy_config &cyl_map_bins(const std::size_t n_phi,
                                       const std::size_t n_z) {
        m_cyl_map_bins = {n_phi, n_z};
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
    constexpr toy_config &beampipe_mat_thickness(const scalar_t t) {
        assert(t > 0.f);
        m_beampipe_mat_thickness = t;
        return *this;
    }
    constexpr toy_config &module_mat_thickness(const scalar_t t) {
        assert(t > 0.f);
        m_module_mat_thickness = t;
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
    constexpr const auto &outer_radius() const { return m_outer_radius; }
    constexpr scalar_t envelope() const { return m_portal_envelope; }
    constexpr scalar_t beampipe_vol_radius() const {
        return m_beampipe_volume_radius;
    }
    constexpr auto &material_config() { return m_material_config; }
    constexpr const auto &material_config() const { return m_material_config; }
    constexpr bool use_material_maps() const { return m_use_material_maps; }
    constexpr const std::array<std::size_t, 2> &cyl_map_bins() const {
        return m_cyl_map_bins;
    }
    constexpr const std::array<std::size_t, 2> &disc_map_bins() const {
        return m_disc_map_bins;
    }
    constexpr scalar_t thickness() const { return m_thickness; }
    constexpr scalar_t beampipe_mat_thickness() const {
        return m_beampipe_mat_thickness;
    }
    constexpr scalar_t module_mat_thickness() const {
        return m_module_mat_thickness;
    }
    auto barrel_mat_generator() const { return m_cyl_mat_generator; }
    auto edc_mat_generator() const { return m_disc_mat_generator; }
    constexpr material<scalar_t> mapped_material() const {
        return m_mapped_material;
    }
    constexpr const auto &barrel_layer_radii() const {
        return m_barrel_layer_radii;
    }
    constexpr const auto &endcap_layer_positions() const {
        return m_endcap_layer_positions;
    }
    constexpr const auto &barrel_layer_binning() const {
        return m_barrel_binning;
    }
    constexpr barrel_generator_config<scalar_t> &barrel_config() {
        return m_barrel_factory_cfg;
    }
    constexpr endcap_generator_config<scalar_t> &endcap_config() {
        return m_endcap_factory_cfg;
    }
    constexpr bool do_check() const { return m_do_check; }
    /// @}
};

namespace detail {

// Helper type, used to define the r- or z-extent of detector volumes
template <typename scalar_t>
struct extent2D {
    scalar_t lower, upper;
};

/// Add the cylinder and disc portals for a volume from explicit parameters
///
/// @param v_builder the volume builder to add the portals to
/// @param names name map for volumes of the detector under construction
/// @param vol_z half length of the cylinder
/// @param h_z half length of the cylinder
/// @param inner_r inner volume radius
/// @param outer_r outer volume radius
/// @param link_north portal volume link of the outer cylinder
/// @param link_south portal volume link of the inner cylinder
/// @param link_east portal volume link of the left disc
/// @param link_west portal volume link of the right disc
template <typename detector_t>
void add_gap_portals(volume_builder_interface<detector_t> *v_builder,
                     typename detector_t::name_map &names,
                     const typename detector_t::scalar_type lower_z,
                     const typename detector_t::scalar_type upper_z,
                     const typename detector_t::scalar_type inner_r,
                     const typename detector_t::scalar_type outer_r,
                     const dindex link_north, const dindex link_south,
                     const dindex link_east, const dindex link_west) {

    using transform3_t = typename detector_t::transform3;
    using scalar_t = typename detector_t::scalar_type;
    using point3_t = typename detector_t::point3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    const transform3_t identity{};
    const dindex vol_idx{v_builder->vol_index()};

    scalar_t min_r{math::min(inner_r, outer_r)};
    scalar_t max_r{math::max(inner_r, outer_r)};
    scalar_t min_z{math::min(lower_z, upper_z)};
    scalar_t max_z{math::max(lower_z, upper_z)};

    using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
    using disc_factory_t = surface_factory<detector_t, ring2D>;

    auto pt_cyl_factory = std::make_shared<cyl_factory_t>();
    auto pt_disc_factory = std::make_shared<disc_factory_t>();

    // Inner cylinder portal
    pt_cyl_factory->push_back({surface_id::e_portal, identity,
                               static_cast<nav_link_t>(link_south),
                               std::vector<scalar_t>{min_r, min_z, max_z}});
    // Outer cylinder portal
    pt_cyl_factory->push_back({surface_id::e_portal, identity,
                               static_cast<nav_link_t>(link_north),
                               std::vector<scalar_t>{max_r, min_z, max_z}});

    // Left disc portal
    pt_disc_factory->push_back({surface_id::e_portal,
                                transform3_t{point3_t{0.f, 0.f, min_z}},
                                static_cast<nav_link_t>(link_west),
                                std::vector<scalar_t>{min_r, max_r}});
    // Right disc portal
    pt_disc_factory->push_back({surface_id::e_portal,
                                transform3_t{point3_t{0.f, 0.f, max_z}},
                                static_cast<nav_link_t>(link_east),
                                std::vector<scalar_t>{min_r, max_r}});

    v_builder->add_surfaces(pt_cyl_factory);
    v_builder->add_surfaces(pt_disc_factory);

    names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);
}

/// Helper method to decorate a volume builder and surface factory with material
///
/// @param cfg config for the toy detector
/// @param det_builder detector builder the volume belongs to
/// @param v_builder the builder of the volume that should be decorated
/// @param sf_factory surface factory that should be decorated with material
///
/// @returns the decorated volume builder and surface factory
template <typename detector_builder_t, typename detector_t>
std::tuple<volume_builder_interface<detector_t> *,
           std::shared_ptr<surface_factory_interface<detector_t>>>
decorate_material(
    const toy_config<typename detector_t::scalar_type> &cfg,
    detector_builder_t &det_builder,
    volume_builder_interface<detector_t> *v_builder,
    std::unique_ptr<
        surface_factory_interface<typename detector_builder_t::detector_type>>
        sf_factory = nullptr) {

    static_assert(
        std::is_same_v<detector_t, typename detector_builder_t::detector_type>,
        "Detector builder and volume builder/surface factory have different "
        "detector type");

    // Decorate the builders with homogeneous material
    if (!cfg.use_material_maps()) {
        // Build the volume with a homogeneous material description
        auto vm_builder =
            det_builder
                .template decorate<homogeneous_material_builder<detector_t>>(
                    v_builder);
        // How to generate the specific material for every surface
        auto mat_generator =
            std::make_shared<homogeneous_material_generator<detector_t>>(
                std::move(sf_factory), cfg.material_config());

        return std::make_tuple(vm_builder, std::move(mat_generator));
    } else {
        // @TODO Add material maps builder here soon
        return std::make_tuple(v_builder, std::move(sf_factory));
    }
}

/// Helper method for creating the barrel surface grids.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param cfg config for the toy detector
/// @param vol_index index of the volume to which the grid should be added
template <typename detector_builder_t>
inline void add_cylinder_grid(
    detector_builder_t &det_builder,
    toy_config<typename detector_builder_t::detector_type::scalar_type> &cfg,
    const dindex vol_index) {

    using detector_t = typename detector_builder_t::detector_type;
    using scalar_t = typename detector_t::scalar_type;

    constexpr auto grid_id = detector_t::accel::id::e_cylinder2_grid;

    using cyl_grid_t =
        typename detector_t::accelerator_container::template get_type<grid_id>;
    using grid_builder_t =
        grid_builder<detector_t, cyl_grid_t, detray::fill_by_pos>;

    const auto &barrel_cfg{cfg.barrel_config()};
    const scalar_t h_z{barrel_cfg.half_length()};

    auto v_builder = det_builder.template decorate<grid_builder_t>(vol_index);
    auto vgr_builder = dynamic_cast<grid_builder_t *>(v_builder);

    vgr_builder->set_type(detector_t::geo_obj_ids::e_sensitive);
    vgr_builder->init_grid(
        {-constant<scalar_t>::pi, constant<scalar_t>::pi, -h_z, h_z},
        {barrel_cfg.binning().first, barrel_cfg.binning().second});
}

/// Helper method for creating the endcap surface grids.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param cfg config for the toy detector
/// @param vol_index index of the volume to which the grid should be added
template <typename detector_builder_t>
inline void add_disc_grid(
    detector_builder_t &det_builder,
    toy_config<typename detector_builder_t::detector_type::scalar_type> &cfg,
    const dindex vol_index) {

    using detector_t = typename detector_builder_t::detector_type;
    using scalar_t = typename detector_t::scalar_type;

    constexpr auto grid_id = detector_t::accel::id::e_disc_grid;

    using disc_grid_t =
        typename detector_t::accelerator_container::template get_type<grid_id>;
    using grid_builder_t =
        grid_builder<detector_t, disc_grid_t, detray::fill_by_pos>;

    const auto &endcap_cfg{cfg.endcap_config()};
    const scalar_t inner_r{cfg.beampipe_vol_radius()};
    const scalar_t outer_r{cfg.outer_radius()};

    auto v_builder = det_builder.template decorate<grid_builder_t>(vol_index);
    auto vgr_builder = dynamic_cast<grid_builder_t *>(v_builder);

    vgr_builder->set_type(detector_t::geo_obj_ids::e_sensitive);
    vgr_builder->init_grid(
        {inner_r, outer_r, -constant<scalar_t>::pi, constant<scalar_t>::pi},
        {endcap_cfg.binning().size(), endcap_cfg.binning().back()});
}

/// Helper method for creating the barrel section.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param gctx geometry context
/// @param cfg config for the toy detector
/// @param names name map for volumes of the detector under construction
/// @param beampipe_idx index of the beampipe outermost volume
///
/// @returns the radial extents of the barrel module layers and gap volumes
template <typename detector_builder_t>
inline auto add_barrel_detector(
    detector_builder_t &det_builder,
    typename detector_builder_t::detector_type::geometry_context &gctx,
    toy_config<typename detector_builder_t::detector_type::scalar_type> &cfg,
    typename detector_builder_t::detector_type::name_map &names,
    dindex beampipe_idx) {

    using detector_t = typename detector_builder_t::detector_type;
    using scalar_t = typename detector_t::scalar_type;
    using transform3_t = typename detector_t::transform3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    // Register the sizes in z per volume index
    std::vector<std::pair<dindex, extent2D<scalar_t>>> volume_sizes{};

    // Mask volume link for portals that exit the detector
    constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
    const transform3_t identity{};

    // Links to the connector gaps. Their volume index depends on the
    // number of endcap layer that are being constructed (before and after
    // the barrel volumes are constructed)
    auto link_east{end_of_world};
    auto link_west{end_of_world};

    if (cfg.n_edc_layers() > 0) {
        link_east = static_cast<nav_link_t>(det_builder.n_volumes() +
                                            2u * cfg.n_brl_layers() + 2u);
        link_west = static_cast<nav_link_t>(det_builder.n_volumes() -
                                            2u * cfg.n_edc_layers() + 1u);
    }

    const scalar_t h_z{cfg.barrel_config().half_length()};
    // Set the inner radius of the first gap to the radius of the beampipe vol.
    scalar_t gap_inner_r{cfg.beampipe_vol_radius()};

    typename cylinder_portal_generator<detector_t>::boundaries vol_bounds{};

    // Alternate barrel module layers and gap volumes
    bool is_gap = true;
    for (unsigned int i = 0u; i < 2u * cfg.n_brl_layers(); ++i) {

        // New volume
        auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
        const dindex vol_idx{v_builder->vol_index()};

        // The barrel volumes are centered at the origin
        v_builder->add_volume_placement(identity);

        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {

            auto link_north{vol_idx - 1};
            auto link_south{vol_idx - 3};

            // The first time a gap is built, it needs to link to the beampipe
            link_south = (i == 1u) ? beampipe_idx : link_south;

            detail::add_gap_portals(v_builder, names, -h_z, h_z, gap_inner_r,
                                    vol_bounds.inner_radius, link_north,
                                    link_south, link_east, link_west);
            volume_sizes.push_back(
                {vol_idx, {gap_inner_r, vol_bounds.inner_radius}});

            // Set the inner gap radius for the next gap volume
            gap_inner_r = vol_bounds.outer_radius;

            names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);

        } else {
            // Limit to maximum valid link
            auto link_north{
                std::min(vol_idx + 3,
                         2 * cfg.n_edc_layers() + 2 * cfg.n_brl_layers() + 1)};
            auto link_south{vol_idx + 1};

            // Configure the module factory for this layer
            auto &barrel_cfg = cfg.barrel_config();

            const unsigned int j{(i + 2u) / 2u};
            barrel_cfg.binning(cfg.barrel_layer_binning().at(j))
                .radius(cfg.barrel_layer_radii().at(j));

            // Configure the portal factory
            cylinder_portal_config<scalar_t> portal_cfg{};

            portal_cfg.envelope(cfg.envelope())
                .fixed_half_length(h_z)
                // Link the volume portals to its neighbors
                .link_north(link_north)
                .link_south(link_south)
                .link_east(link_east)
                .link_west(link_west);

            // Configure the material
            cfg.material_config().thickness(cfg.module_mat_thickness());

            // Add a layer of module surfaces
            auto module_factory =
                std::make_unique<barrel_generator<detector_t, rectangle2D>>(
                    barrel_cfg);

            // Add cylinder and disc portals
            auto portal_factory =
                std::make_shared<cylinder_portal_generator<detector_t>>(
                    portal_cfg);

            auto [vm_builder, module_mat_factory] = decorate_material(
                cfg, det_builder, v_builder, std::move(module_factory));

            vm_builder->add_surfaces(module_mat_factory, gctx);
            vm_builder->add_surfaces(portal_factory);

            // Set the new current boundaries, to construct the next gap
            vol_bounds = portal_factory->volume_boundaries();
            volume_sizes.push_back(
                {vol_idx, {vol_bounds.inner_radius, vol_bounds.outer_radius}});

            names[vol_idx + 1u] = "barrel_" + std::to_string(vol_idx);

            // Add a cylinder grid to every barrel module layer
            add_cylinder_grid(det_builder, cfg, vol_idx);
        }
    }

    // Add a final gap volume to get to the full barrel radius
    auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
    const dindex vol_idx{v_builder->vol_index()};
    v_builder->add_volume_placement(identity);

    detail::add_gap_portals(v_builder, names, -h_z, h_z,
                            vol_bounds.outer_radius, cfg.outer_radius(),
                            end_of_world, vol_idx - 2u, link_east, link_west);
    volume_sizes.push_back(
        {vol_idx, {vol_bounds.outer_radius, cfg.outer_radius()}});

    names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);

    return volume_sizes;
}

/// Helper method for creating one of the two endcaps.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param gctx geometry context
/// @param cfg config for the toy detector
/// @param names name map for volumes of the detector under construction
/// @param beampipe_idx index of the beampipe outermost volume
///
/// @returns the z extents of the endcap module layers and gap volumes
template <typename detector_builder_t>
inline auto add_endcap_detector(
    detector_builder_t &det_builder,
    typename detector_builder_t::detector_type::geometry_context &gctx,
    toy_config<typename detector_builder_t::detector_type::scalar_type> &cfg,
    typename detector_builder_t::detector_type::name_map &names,
    dindex beampipe_idx) {

    using detector_t = typename detector_builder_t::detector_type;
    using scalar_t = typename detector_t::scalar_type;
    using point3_t = typename detector_t::point3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    std::vector<std::pair<dindex, extent2D<scalar_t>>> volume_sizes{};

    // Mask volume link for portals that exit the detector
    constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
    // Disc portal link to the barrel section (via connector gap volume, which
    // is the second volume initialized by this function, but built in a
    // later step, when all necessaty information is available)
    const dindex connector_link{det_builder.n_volumes() + 1u};

    const auto sign{static_cast<scalar_t>(cfg.endcap_config().side())};

    // Mask volume links
    auto link_north{end_of_world};
    auto link_south{beampipe_idx};

    // Set the inner radius of the first gap to the radius of the beampipe vol.
    const scalar_t inner_radius{cfg.beampipe_vol_radius()};
    const scalar_t outer_radius{cfg.outer_radius()};

    scalar_t gap_east_z{sign * cfg.barrel_config().half_length()};

    typename cylinder_portal_generator<detector_t>::boundaries vol_bounds{};

    // Alternate endcap module layers and gap volumes
    bool is_gap = true;
    for (dindex i = 0u; i < 2u * cfg.n_edc_layers(); ++i) {

        // New volume
        auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
        const dindex vol_idx{v_builder->vol_index()};

        // Don't build the first gap here, as this will be done in a separate
        // step once all portals of the barrel are constructed
        if (i == 1u) {
            const scalar_t gap_west_z{sign *
                                      std::min(std::abs(vol_bounds.upper_z),
                                               std::abs(vol_bounds.lower_z))};

            volume_sizes.push_back({vol_idx, {gap_east_z, gap_west_z}});

            // Update the gap extent for the next gap
            gap_east_z = sign * std::max(std::abs(vol_bounds.upper_z),
                                         std::abs(vol_bounds.lower_z));

            names[vol_idx + 1u] = "connector_gap_" + std::to_string(vol_idx);

            is_gap = !is_gap;
            continue;
        }

        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {
            auto link_east{static_cast<dindex>(static_cast<int>(vol_idx) +
                                               cfg.endcap_config().side() - 2)};
            auto link_west{static_cast<dindex>(static_cast<int>(vol_idx) -
                                               cfg.endcap_config().side() - 2)};

            const scalar_t gap_west_z{sign *
                                      std::min(std::abs(vol_bounds.upper_z),
                                               std::abs(vol_bounds.lower_z))};

            const point3_t gap_center{0.f, 0.f,
                                      0.5f * (gap_east_z + gap_west_z)};
            v_builder->add_volume_placement({gap_center});

            detail::add_gap_portals(v_builder, names, gap_west_z, gap_east_z,
                                    inner_radius, outer_radius, link_north,
                                    link_south, link_east, link_west);

            volume_sizes.push_back({vol_idx, {gap_east_z, gap_west_z}});

            // Set the inner gap radius for the next gap volume
            gap_east_z = sign * std::max(std::abs(vol_bounds.upper_z),
                                         std::abs(vol_bounds.lower_z));

            names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);

        } else {
            const dindex j{i / 2u};

            auto link_east{static_cast<dindex>(static_cast<int>(vol_idx) +
                                               cfg.endcap_config().side() + 2)};
            auto link_west{static_cast<dindex>(static_cast<int>(vol_idx) -
                                               cfg.endcap_config().side() + 2)};

            // The first endacp layer needs to link to the connector gap
            // The last endcap layer needs to exit the detector
            if (sign < 0) {
                link_east = (i == 0u) ? connector_link : link_east;
                link_west =
                    (j == cfg.n_edc_layers() - 1u) ? end_of_world : link_west;
            } else {
                link_west = (i == 0u) ? connector_link : link_west;
                link_east =
                    (j == cfg.n_edc_layers() - 1u) ? end_of_world : link_east;
            }

            // Position the volume at the respective endcap layer position
            const scalar_t center_z{sign * cfg.endcap_layer_positions().at(j)};
            const point3_t vol_center{0.f, 0.f, center_z};
            v_builder->add_volume_placement({vol_center});

            // Configure the module factory for this layer
            auto &endcap_cfg = cfg.endcap_config();

            endcap_cfg.center(center_z)
                .inner_radius(inner_radius)
                .outer_radius(outer_radius);

            // Configure the portal factory
            cylinder_portal_config<scalar_t> portal_cfg{};

            portal_cfg.envelope(cfg.envelope())
                .fixed_inner_radius(inner_radius)
                .fixed_outer_radius(outer_radius)
                // Link the volume portals to their neighbors
                .link_north(link_north)
                .link_south(link_south)
                .link_east(link_east)
                .link_west(link_west);

            // Configure the material
            cfg.material_config().thickness(cfg.module_mat_thickness());

            // Add a layer of module surfaces
            auto module_factory =
                std::make_unique<endcap_generator<detector_t, trapezoid2D>>(
                    endcap_cfg);

            // Add cylinder and disc portals
            auto portal_factory =
                std::make_shared<cylinder_portal_generator<detector_t>>(
                    portal_cfg);

            auto [vm_builder, module_mat_factory] = decorate_material(
                cfg, det_builder, v_builder, std::move(module_factory));

            vm_builder->add_surfaces(module_mat_factory, gctx);
            vm_builder->add_surfaces(portal_factory);

            // Set the new current boundaries to construct the next gap
            vol_bounds = portal_factory->volume_boundaries();
            volume_sizes.push_back(
                {vol_idx, {vol_bounds.lower_z, vol_bounds.upper_z}});

            names[vol_idx + 1u] = "endcap_" + std::to_string(vol_idx);

            // Add a disc grid to every endcap module layer
            add_disc_grid(det_builder, cfg, vol_idx);
        }
    }
    return volume_sizes;
}

/// Helper method for creating the portals of one of endcap sections of the
/// beampipe volume.
///
/// @param beampipe_builder volume builder for the beampipe volume
/// @param cfg config for the toy detector
/// @param neg_edc_lay_sizes indices and z-extent of the endcap volumes of one
///                          detector side (positive or negative)
template <typename detector_builder_t, typename vol_extent_data_t>
inline void add_connector_portals(
    detector_builder_t &det_builder,
    toy_config<typename detector_builder_t::detector_type::scalar_type> &cfg,
    const dindex beampipe_idx, const vol_extent_data_t edc_vol_extents,
    const vol_extent_data_t &brl_vol_extents) {

    using detector_t = typename detector_builder_t::detector_type;
    using transform3_t = typename detector_t::transform3;
    using scalar_t = typename detector_t::scalar_type;
    using point3_t = typename detector_t::point3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
    using disc_factory_t = surface_factory<detector_t, ring2D>;

    auto pt_cyl_factory = std::make_shared<cyl_factory_t>();
    auto pt_disc_factory = std::make_shared<disc_factory_t>();

    // Mask volume link for portals that exit the detector
    constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
    const transform3_t identity{};

    // Set the inner radius of the gap to the radius of the beampipe vol.
    const scalar_t inner_r{cfg.beampipe_vol_radius()};
    const scalar_t outer_r{cfg.outer_radius()};

    // Get the data of the gap volume
    // The connector gap is always the second volume constructed in the endcaps
    const auto &connector_gap_data = edc_vol_extents[1];
    const dindex connector_gap_idx{connector_gap_data.first};

    const scalar_t side{std::copysign(1.f, connector_gap_data.second.lower)};
    const scalar_t gap_east_z{connector_gap_data.second.lower};
    const scalar_t gap_west_z{connector_gap_data.second.upper};
    const scalar_t min_z{math::min(gap_east_z, gap_west_z)};
    const scalar_t max_z{math::max(gap_east_z, gap_west_z)};
    const point3_t gap_center{0.f, 0.f, 0.5f * (gap_east_z + gap_west_z)};

    // Check that the volume under construction is really the connector gap
    assert(std::abs(gap_east_z) == cfg.barrel_config().half_length() ||
           std::abs(gap_west_z) == cfg.barrel_config().half_length());

    volume_builder_interface<detector_t> *connector_builder =
        det_builder[connector_gap_idx];

    connector_builder->add_volume_placement({gap_center});

    // Go over all volume extends and build the corresponding disc portal
    for (const auto &e : brl_vol_extents) {

        const scalar_t min_r{math::min(e.second.lower, e.second.upper)};
        const scalar_t max_r{math::max(e.second.lower, e.second.upper)};

        // Barrel-facing disc portal
        pt_disc_factory->push_back(
            {surface_id::e_portal,
             transform3_t{point3_t{0.f, 0.f, (side < 0.f) ? max_z : min_z}},
             static_cast<nav_link_t>(e.first),
             std::vector<scalar_t>{min_r, max_r}});
    }

    // Outward-facing disc portal
    pt_disc_factory->push_back(
        {surface_id::e_portal,
         transform3_t{point3_t{0.f, 0.f, (side < 0.f) ? min_z : max_z}},
         static_cast<nav_link_t>(connector_gap_idx - 1u),
         std::vector<scalar_t>{inner_r, outer_r}});

    // Inner cylinder portal
    pt_cyl_factory->push_back({surface_id::e_portal, identity,
                               static_cast<nav_link_t>(beampipe_idx),
                               std::vector<scalar_t>{inner_r, min_z, max_z}});

    // Outer cylinder portal
    pt_cyl_factory->push_back({surface_id::e_portal, identity, end_of_world,
                               std::vector<scalar_t>{outer_r, min_z, max_z}});

    // Add all portals to the beampipe volume
    connector_builder->add_surfaces(pt_disc_factory);
    connector_builder->add_surfaces(pt_cyl_factory);
}

/// Helper method for creating the portals of the beampipe barrel section.
///
/// @param beampipe_builder volume builder for the beampipe volume
/// @param cfg config for the toy detector
template <typename detector_t>
inline void add_beampipe_portals(
    volume_builder_interface<detector_t> *beampipe_builder,
    toy_config<typename detector_t::scalar_type> &cfg) {

    using scalar_t = typename detector_t::scalar_type;
    using transform3_t = typename detector_t::transform3;
    using point3_t = typename detector_t::point3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
    using disc_factory_t = surface_factory<detector_t, ring2D>;

    auto pt_cyl_factory = std::make_shared<cyl_factory_t>();

    // Mask volume link for portals that exit the detector
    constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};

    const scalar_t h_z{cfg.barrel_config().half_length()};
    const scalar_t inner_r{0.f};
    const scalar_t outer_r{cfg.beampipe_vol_radius()};

    // If there are no endcaps, add the beampipe disc portals at the boundaries
    // of the barrel section
    if (cfg.n_edc_layers() == 0u) {

        auto pt_disc_factory = std::make_shared<disc_factory_t>();

        // Lower dics portal
        pt_disc_factory->push_back(
            {surface_id::e_portal, transform3_t{point3_t{0.f, 0.f, -h_z}},
             end_of_world, std::vector<scalar_t>{inner_r, outer_r}});

        // Outward-facing disc portal
        pt_disc_factory->push_back(
            {surface_id::e_portal, transform3_t{point3_t{0.f, 0.f, h_z}},
             end_of_world, std::vector<scalar_t>{inner_r, outer_r}});

        beampipe_builder->add_surfaces(pt_disc_factory);
    }

    // Cylinder portal that leads into the barrel section
    dindex first_barrel_idx{
        cfg.n_edc_layers() == 0u ? end_of_world : 2u * cfg.n_edc_layers() + 2u};
    pt_cyl_factory->push_back(
        {surface_id::e_portal, transform3_t{},
         static_cast<nav_link_t>(first_barrel_idx),
         std::vector<scalar_t>{
             outer_r, -h_z + std::numeric_limits<scalar_t>::epsilon(),
             h_z - std::numeric_limits<scalar_t>::epsilon()}});

    beampipe_builder->add_surfaces(pt_cyl_factory);
}

/// Helper method for creating the portals of one of endcap sections of the
/// beampipe volume.
///
/// @param beampipe_builder volume builder for the beampipe volume
/// @param cfg config for the toy detector
/// @param neg_edc_lay_sizes indices and z-extent of the endcap volumes of one
///                          detector side (positive or negative)
template <typename detector_t, typename scalar_t, typename layer_size_cont_t>
inline void add_beampipe_portals(
    volume_builder_interface<detector_t> *beampipe_builder,
    toy_config<scalar_t> &cfg, const layer_size_cont_t &edc_lay_sizes) {

    using transform3_t = typename detector_t::transform3;
    using point3_t = typename detector_t::point3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;

    using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
    using disc_factory_t = surface_factory<detector_t, ring2D>;

    // Mask volume link for portals that exit the detector
    constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
    // For which endcap side should the portals be constructed?
    const scalar_t side{std::copysign(1.f, edc_lay_sizes.front().second.lower)};

    auto pt_cyl_factory = std::make_shared<cyl_factory_t>();
    auto pt_disc_factory = std::make_shared<disc_factory_t>();

    const scalar_t inner_r{0.f};
    const scalar_t outer_r{cfg.beampipe_vol_radius()};
    scalar_t disc_pos_z{-side * std::numeric_limits<scalar_t>::max()};

    // Go over all volume extends and build the corresponding cylinder portal
    for (const auto &e : edc_lay_sizes) {

        const scalar_t min_z{math::min(e.second.lower, e.second.upper)};
        const scalar_t max_z{math::max(e.second.lower, e.second.upper)};

        if (side < 0.f) {
            disc_pos_z = (disc_pos_z > min_z) ? min_z : disc_pos_z;
        } else {
            disc_pos_z = (disc_pos_z < max_z) ? max_z : disc_pos_z;
        }

        pt_cyl_factory->push_back(
            {surface_id::e_portal, transform3_t{},
             static_cast<nav_link_t>(e.first),
             std::vector<scalar_t>{outer_r, min_z, max_z}});
    }

    // Outward-facing disc portal
    pt_disc_factory->push_back(
        {surface_id::e_portal, transform3_t{point3_t{0.f, 0.f, disc_pos_z}},
         end_of_world, std::vector<scalar_t>{inner_r, outer_r}});

    // Add all portals to the beampipe volume
    beampipe_builder->add_surfaces(pt_disc_factory);
    beampipe_builder->add_surfaces(pt_cyl_factory);
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
    using transform3_t = typename detector_t::transform3;
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
    using vol_extent_container_t =
        std::vector<std::pair<dindex, detail::extent2D<scalar_t>>>;

    static_assert(std::is_same_v<typename detector_t::scalar_type, scalar_t>,
                  "Scalar type used for toy detector config does not match the "
                  "detector algebra type");

    // Check config
    if (cfg.n_edc_layers() > cfg.endcap_layer_positions().size()) {
        throw std::invalid_argument(
            "ERROR: Too many endcap layers requested (max " +
            std::to_string(cfg.endcap_layer_positions().size()) + ")!");
    }
    if (cfg.n_brl_layers() > cfg.barrel_layer_radii().size() - 1u) {
        throw std::invalid_argument(
            "ERROR: Too many barrel layers requested (max " +
            std::to_string(cfg.barrel_layer_radii().size() - 1u) + ")!");
    }
    if (cfg.n_edc_layers() > 0 && cfg.n_brl_layers() < 4) {
        throw std::invalid_argument(
            "ERROR: All four barrel layers need to be present in order to add "
            "endcap layers");
    }

    // Toy detector builder
    builder_t det_builder;

    // Detector and volume names
    typename detector_t::name_map name_map = {{0u, "toy_detector_new"}};
    // Geometry context object
    typename detector_t::geometry_context gctx{};

    // Add the volume that contains the beampipe
    cfg.material_config().thickness(cfg.beampipe_mat_thickness());
    auto [beampipe_builder, pt_cyl_factory] = detail::decorate_material(
        cfg, det_builder, det_builder.new_volume(volume_id::e_cylinder),
        std::make_unique<cyl_factory_t>());

    const dindex beampipe_idx{beampipe_builder->vol_index()};
    beampipe_builder->add_volume_placement(transform3_t{});
    name_map[beampipe_idx + 1u] = "beampipe_" + std::to_string(beampipe_idx);

    // Add the beampipe as a passive material surface
    scalar_t max_z{cfg.n_edc_layers() == 0u ? cfg.barrel_config().half_length()
                                            : cfg.endcap_layer_positions().at(
                                                  cfg.n_edc_layers() - 1u)};
    scalar_t min_z{-max_z};

    pt_cyl_factory->push_back(
        {surface_id::e_passive, transform3_t{},
         static_cast<nav_link_t>(beampipe_idx),
         std::vector<scalar_t>{cfg.barrel_layer_radii().at(0), min_z, max_z}});

    beampipe_builder->add_surfaces(pt_cyl_factory);

    // Build the negative endcap
    vol_extent_container_t neg_edc_vol_extents;
    if (cfg.n_edc_layers() > 0u) {
        cfg.endcap_config().side(-1);

        neg_edc_vol_extents = detail::add_endcap_detector(
            det_builder, gctx, cfg, name_map, beampipe_idx);

        // Add the beampipe volume portals for the negative endcap section
        detail::add_beampipe_portals(beampipe_builder, cfg,
                                     neg_edc_vol_extents);
    }
    // Build the barrel section
    vol_extent_container_t brl_vol_extents;
    if (cfg.n_brl_layers() > 0u) {

        brl_vol_extents = detail::add_barrel_detector(det_builder, gctx, cfg,
                                                      name_map, beampipe_idx);

        // Add the beampipe volume portals for the barrel section
        detail::add_beampipe_portals(beampipe_builder, cfg);
    }
    // Build the positive endcap
    vol_extent_container_t pos_edc_vol_extents;
    if (cfg.n_edc_layers() > 0u) {
        cfg.endcap_config().side(1);

        pos_edc_vol_extents = detail::add_endcap_detector(
            det_builder, gctx, cfg, name_map, beampipe_idx);

        // Add the beampipe volume portals for the positive endcap section
        detail::add_beampipe_portals(beampipe_builder, cfg,
                                     pos_edc_vol_extents);

        // Add the connection between barrel and both endcaps
        // Negative endcap
        add_connector_portals(det_builder, cfg, beampipe_idx,
                              neg_edc_vol_extents, brl_vol_extents);
        // Positive endcap
        add_connector_portals(det_builder, cfg, beampipe_idx,
                              pos_edc_vol_extents, brl_vol_extents);
    }

    // Build and return the detector
    auto det = det_builder.build(resource);

    if (cfg.do_check()) {
        detray::detail::check_consistency(det);
    }

    return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
