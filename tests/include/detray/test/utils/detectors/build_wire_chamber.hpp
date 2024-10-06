/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
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
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/utils/consistency_checker.hpp"

// Detray test include(s)
#include "detray/test/utils/detectors/factories/wire_layer_generator.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <format>

namespace detray {

/// Configuration for building a wire chamber detector
struct wire_chamber_config {

    /// Number of layers
    unsigned int m_n_layers{10u};
    /// The inner radius of the first layer
    scalar m_first_layer_inner_rad{500.f * unit<scalar>::mm};
    /// Half z of cylinder chamber
    scalar m_half_z{1000.f * unit<scalar>::mm};
    /// Radius of the material rods
    scalar m_mat_radius{15.f * unit<scalar>::um};
    /// Type of material for the material rods
    material<scalar> m_wire_mat{tungsten<scalar>()};
    /// Config for the wire generation (barrel)
    wire_layer_generator_config<scalar> m_wire_factory_cfg{};
    /// Do a full detector consistency check after building
    bool m_do_check{true};

    /// Setters
    /// @{
    constexpr wire_chamber_config &n_layers(const unsigned int n) {
        m_n_layers = n;
        return *this;
    }
    constexpr wire_chamber_config &first_layer_inner_radius(const scalar r) {
        m_first_layer_inner_rad = r;
        return *this;
    }
    constexpr wire_chamber_config &half_z(const scalar hz) {
        m_half_z = hz;
        return *this;
    }
    constexpr wire_chamber_config &cell_size(const scalar c) {
        m_wire_factory_cfg.cell_size(c);
        return *this;
    }
    constexpr wire_chamber_config &stereo_angle(const scalar s) {
        m_wire_factory_cfg.stereo_angle(s);
        return *this;
    }
    constexpr wire_chamber_config &mat_radius(const scalar r) {
        m_mat_radius = r;
        return *this;
    }
    constexpr wire_chamber_config &wire_material(const material<scalar> &m) {
        m_wire_mat = m;
        return *this;
    }
    constexpr wire_chamber_config &do_check(const bool check) {
        m_do_check = check;
        return *this;
    }
    /// @}

    /// Getters
    /// @{
    constexpr unsigned int n_layers() const { return m_n_layers; }
    constexpr scalar first_layer_inner_radius() const {
        return m_first_layer_inner_rad;
    }
    constexpr scalar half_z() const { return m_half_z; }
    constexpr scalar cell_size() const {
        return m_wire_factory_cfg.cell_size();
    }
    constexpr scalar stereo_angle() const {
        return m_wire_factory_cfg.stereo_angle();
    }
    constexpr scalar mat_radius() const { return m_mat_radius; }
    constexpr const material<scalar> &wire_material() const {
        return m_wire_mat;
    }
    constexpr wire_layer_generator_config<scalar> &layer_config() {
        return m_wire_factory_cfg;
    }
    constexpr bool do_check() const { return m_do_check; }
    /// @}

    private:
    /// Print the wire chamber configuration
    friend inline std::ostream &operator<<(std::ostream &out,
                                           const wire_chamber_config &cfg) {
        out << "\nWire Chamber\n"
            << "----------------------------\n"
            << "  No. layers            : " << cfg.n_layers() << "\n"
            << "  First layer inner rad.: " << cfg.first_layer_inner_radius()
            << " [mm]\n"
            << "  Half length z         : " << cfg.half_z() << " [mm]\n"
            << "  Cell size             : " << cfg.cell_size() << " [mm]\n"
            << "  Stereo angle          : " << cfg.stereo_angle() << " [rad]\n"
            << "  Wire material         : " << cfg.wire_material() << "\n"
            << "  Material rad.         : " << cfg.mat_radius() << " [mm]\n";

        return out;
    }

};  // wire chamber config

inline auto build_wire_chamber(vecmem::memory_resource &resource,
                               wire_chamber_config &cfg) {

    using builder_t = detector_builder<default_metadata, volume_builder>;
    using detector_t = typename builder_t::detector_type;
    using scalar_t = typename detector_t::scalar_type;

    // Wire chamber detector builder
    builder_t det_builder;

    // Detector and volume names
    typename detector_t::name_map name_map = {{0u, "wire_chamber"}};
    // Geometry context object
    typename detector_t::geometry_context gctx{};

    // Navigation link when leaving the detector
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};
    const scalar_t inner_rad{cfg.first_layer_inner_radius()};
    const scalar_t cell_size{cfg.cell_size()};

    //
    // Build empty inner volume, where silicon subdetectors would sit
    //
    auto inner_v_builder = det_builder.new_volume(volume_id::e_cylinder);
    inner_v_builder->add_volume_placement();
    const dindex inner_vol_idx{inner_v_builder->vol_index()};
    name_map[1u] = "inner_vol_0";

    // Configure the portal factory
    // TODO: Add material maps that model the silicon detector budget
    cylinder_portal_config<scalar_t> inner_pt_cfg{};

    inner_pt_cfg.do_autofit(false)
        .fixed_half_length(cfg.half_z())
        .fixed_inner_radius(0.f)
        .fixed_outer_radius(inner_rad)
        // No inner subdetectors present -> don't build inner portal
        .build_inner(false)
        .link_north(inner_vol_idx + 1u)
        .link_south(leaving_world)
        .link_east(leaving_world)
        .link_west(leaving_world);

    auto inner_pt_factory =
        std::make_shared<cylinder_portal_generator<detector_t>>(inner_pt_cfg);
    inner_v_builder->add_surfaces(inner_pt_factory, gctx);

    //
    // Build layer volumes
    //
    const unsigned int n_layers{cfg.n_layers()};
    for (unsigned int i_lay = 0; i_lay < n_layers; i_lay++) {

        // New volume
        auto vm_builder = det_builder.new_volume(volume_id::e_cylinder);
        // auto vm_builder = decorate_material(cfg, det_builder, v_builder);
        const dindex vol_idx{vm_builder->vol_index()};

        // The barrel volumes are centered at the origin
        vm_builder->add_volume_placement();

        // The maximal inner and outer radius of the volume
        const scalar_t inner_layer_rad =
            inner_rad + static_cast<scalar_t>(i_lay) * 2.f * cell_size;
        const scalar_t outer_layer_rad =
            inner_rad + static_cast<scalar_t>(i_lay + 1u) * 2.f * cell_size;

        // Configure the wire layer factory for this layer
        auto &layer_cfg = cfg.layer_config();
        layer_cfg.inner_layer_radius(inner_layer_rad).half_length(cfg.half_z());

        const scalar_t sign = (i_lay % 2 == 0) ? 1 : -1;
        layer_cfg.stereo_angle(sign * math::fabs(layer_cfg.stereo_angle()));

        // Configure the portal factory
        cylinder_portal_config<scalar_t> layer_portal_cfg{};

        // Limit to maximum valid link
        auto link_north{i_lay == cfg.n_layers() - 1u ? leaving_world
                                                     : vol_idx + 1u};
        auto link_south{vol_idx - 1u};

        layer_portal_cfg.do_autofit(false)
            .fixed_half_length(cfg.half_z())
            .fixed_inner_radius(inner_layer_rad)
            .fixed_outer_radius(outer_layer_rad)
            // Link the volume portals to its neighbors
            .link_north(link_north)
            .link_south(link_south)
            .link_east(leaving_world)
            .link_west(leaving_world);

        // Configure the material
        // cfg.material_config().thickness(cfg.mat_radius());
        auto wire_mat_factory =
            std::make_shared<wire_layer_generator<detector_t>>(layer_cfg);

        auto portal_mat_factory =
            std::make_shared<cylinder_portal_generator<detector_t>>(
                layer_portal_cfg);

        // Add a layer of module surfaces (may have material)
        /*auto wire_mat_factory = decorate_material<detector_t>(
            cfg,
            std::make_unique<wire_layer_generator<detector_t, line_square>>(
                layer_cfg),
            true);

        // Add cylinder and disc portals (may have material, depending on
        // whether material maps are being used or not)
        auto portal_mat_factory = decorate_material<detector_t>(
            cfg, std::make_unique<cylinder_portal_generator<detector_t>>(
                        layer_portal_cfg));*/

        vm_builder->add_surfaces(portal_mat_factory);
        vm_builder->add_surfaces(wire_mat_factory, gctx);

        name_map[vol_idx + 1u] = std::format("layer_vol_{}", vol_idx);

        // Add a cylinder grid to every barrel module layer
        // add_cylinder_grid(det_builder, cfg, vol_idx);
    }

    // Build and return the detector
    auto det = det_builder.build(resource);

    if (cfg.do_check()) {
        const bool verbose_check{false};
        detray::detail::check_consistency(det, verbose_check, name_map);
    }

    return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
