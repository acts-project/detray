/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/telescope_metadata.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tools/cuboid_portal_generator.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/material_builder.hpp"
#include "detray/tools/telescope_generator.hpp"
#include "detray/utils/consistency_checker.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

namespace detray {

namespace {

/// Configure the toy detector
template <typename mask_shape_t = rectangle2D<>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>>
struct tel_det_config {

    using vector3 = __plugin::vector3<detray::scalar>;

    /// Construct from existing mask
    tel_det_config(const mask<mask_shape_t> &m, const trajectory_t &t = {})
        : m_mask(m), m_trajectory(t) {}

    /// Construct from mask parameters (except volume link, which is not needed)
    template <
        typename... Args,
        std::enable_if_t<(std::is_same_v<Args, scalar> || ...), bool> = true>
    tel_det_config(Args &&... args) : m_mask(0u, std::forward<Args>(args)...) {}

    /// Mask of the test surfaces
    mask<mask_shape_t> m_mask;
    /// No. of test surfaces in the telescope
    unsigned int m_n_surfaces{10u};
    /// Length of the telescope
    scalar m_length{500.f * unit<scalar>::mm};
    /// Concrete positions where to place the surfaces along the pilot track
    std::vector<scalar> m_positions{};
    /// Material for the test surfaces
    material<scalar> m_material = silicon_tml<scalar>();
    /// Thickness of the material
    scalar m_thickness{80.f * unit<scalar>::um};
    /// Pilot track along which to place the surfaces
    trajectory_t m_trajectory{};
    /// Safety envelope between the test surfaces and the portals
    scalar m_envelope{0.1f * unit<scalar>::mm};
    /// Run detector consistency check after reading
    bool m_do_check{true};

    /// Setters
    /// @{
    constexpr tel_det_config &module_mask(const mask<mask_shape_t> &m) {
        m_mask = m;
        return *this;
    }
    constexpr tel_det_config &n_surfaces(const unsigned int n) {
        m_n_surfaces = n;
        return *this;
    }
    constexpr tel_det_config &length(const scalar l) {
        assert((l > 0.f) &&
               "Telescope detector length must be greater than zero");
        m_length = l;
        return *this;
    }
    constexpr tel_det_config &positions(const std::vector<scalar> &dists) {
        m_positions.clear();
        std::copy_if(dists.begin(), dists.end(),
                     std::back_inserter(m_positions),
                     [](scalar d) { return (d >= 0.f); });
        return *this;
    }
    constexpr tel_det_config &module_material(const material<scalar> &mat) {
        m_material = mat;
        return *this;
    }
    constexpr tel_det_config &mat_thickness(const scalar t) {
        assert(t > 0.f && "Material thickness must be greater than zero");
        m_thickness = t;
        return *this;
    }
    constexpr tel_det_config &pilot_track(const trajectory_t &traj) {
        m_trajectory = traj;
        return *this;
    }
    constexpr tel_det_config &envelope(const scalar e) {
        assert(e > 0.f && "Portal envelope must be greater than zero");
        m_envelope = e;
        return *this;
    }
    tel_det_config &do_check(const bool check) {
        m_do_check = check;
        return *this;
    }
    /// @}

    /// Getters
    /// @{
    constexpr const mask<mask_shape_t> &module_mask() const { return m_mask; }
    constexpr unsigned int n_surfaces() const { return m_n_surfaces; }
    constexpr scalar length() const { return m_length; }
    const std::vector<scalar> &positions() const { return m_positions; }
    constexpr const material<scalar> &module_material() const {
        return m_material;
    }
    constexpr scalar mat_thickness() const { return m_thickness; }
    const trajectory_t &pilot_track() const { return m_trajectory; }
    constexpr scalar envelope() const { return m_envelope; }
    bool do_check() const { return m_do_check; }
    /// @}
};

}  // namespace

/// Builds a detray geometry that contains only one volume with one type of
/// surfaces, where the last surface is the portal that leaves the telescope.
/// The detector is auto-constructed by following a trajectory state through
/// space. The trajectory and the given distances determine the positions of
/// the plane surfaces. The dis
///
/// @tparam mask_t the type of mask for the telescope surfaces
/// @tparam trajectory_t the type of the pilot trajectory
///
/// @param resource the memory resource for the detector containers
/// @param cfg configuration struct of the telescope detector
///
/// @returns a complete detector object
template <typename mask_shape_t = rectangle2D<>,
          typename trajectory_t = detail::ray<__plugin::transform3<scalar>>>
inline auto create_telescope_detector(
    vecmem::memory_resource &resource,
    const tel_det_config<mask_shape_t, trajectory_t> &cfg = {
        20.f * unit<scalar>::mm, 20.f * unit<scalar>::mm}) {

    detector_builder<telescope_metadata<mask_shape_t>, volume_builder>
        det_builder;

    using detector_t = typename decltype(det_builder)::detector_type;

    // Detector and volume names
    typename detector_t::name_map name_map = {{0u, "telescope_detector"},
                                              {1u, "telescope_world_0"}};

    // Create an empty cuboid volume with homogeneous material description
    auto v_builder = det_builder.new_volume(volume_id::e_cuboid);
    const dindex vol_idx{v_builder->vol_index()};
    auto vm_builder =
        det_builder.template decorate<material_builder<detector_t>>(vol_idx);

    // Identity transform
    vm_builder->add_volume_placement();

    // Add module surfaces to volume
    using telescope_factory =
        telescope_generator<detector_t, mask_shape_t, trajectory_t>;
    auto tel_generator =
        cfg.positions().empty()
            ? std::make_unique<telescope_factory>(
                  cfg.length(), cfg.n_surfaces(), cfg.module_mask().values(),
                  cfg.pilot_track())
            : std::make_unique<telescope_factory>(cfg.positions(),
                                                  cfg.module_mask().values(),
                                                  cfg.pilot_track());

    // @todo: Temporal restriction due to missing local navigation
    assert((tel_generator->size() < 20u) &&
           "Due to WIP, please choose less than 20 surfaces for now");

    std::vector<material_data<scalar>> sf_materials(
        tel_generator->size(),
        material_data<scalar>{cfg.mat_thickness(), cfg.module_material()});

    using material_id = typename detector_t::materials::id;
    constexpr bool is_line{std::is_same_v<mask_shape_t, detray::line<true>> ||
                           std::is_same_v<mask_shape_t, detray::line<false>>};
    const auto mat_id = is_line ? material_id::e_rod : material_id::e_slab;

    auto tel_mat_generator = std::make_shared<material_factory<detector_t>>(
        std::move(tel_generator));
    tel_mat_generator->add_material(mat_id, std::move(sf_materials));

    // Add a portal box around the cuboid volume
    auto portal_generator = std::make_shared<material_factory<detector_t>>(
        std::make_unique<cuboid_portal_generator<detector_t>>(cfg.envelope()));

    // @TODO: Put no material instead of 'vacuum'
    std::vector<material_data<scalar>> pt_materials(
        portal_generator->size(), material_data<scalar>{0.f, vacuum<scalar>{}});
    portal_generator->add_material(detector_t::materials::id::e_slab,
                                   std::move(pt_materials));

    vm_builder->add_surfaces(tel_mat_generator);
    vm_builder->add_surfaces(portal_generator);

    det_builder.set_volume_finder(resource);
    det_builder.volume_finder().push_back(std::vector<dindex>{vol_idx});

    // Build and return the detector
    auto det = det_builder.build(resource);

    if (cfg.do_check()) {
        detray::detail::check_consistency(det);
    }

    return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
