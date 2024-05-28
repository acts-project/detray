/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/utils/create_path.hpp"
#include "detray/io/utils/file_handle.hpp"
#include "detray/materials/detail/material_accessor.hpp"
#include "detray/propagator/base_actor.hpp"

// System include(s)
#include <filesystem>

namespace detray::material_validator {

/// @brief Record the material budget per thickness or pathlength
template <typename scalar_t>
struct material_record {
    /// Phi and eta values of the track for which the material was recorded
    /// @{
    scalar_t phi{detail::invalid_value<scalar_t>()};
    scalar_t eta{detail::invalid_value<scalar_t>()};
    /// @}
    /// Accumulated radiation length per pathlength through the material
    scalar_t sX0{0.f};
    /// Accumulated radiation length per thickness
    scalar_t tX0{0.f};
    /// Accumulated interaction length per pathlength through the material
    scalar_t sL0{0.f};
    /// Accumulated interaction length per thickness
    scalar_t tL0{0.f};
};

/// @brief Functor to retrieve the material parameters for a given local
/// position
struct get_material_params {

    template <typename mat_group_t, typename index_t, typename point2_t,
              typename scalar_t>
    inline auto operator()(
        [[maybe_unused]] const mat_group_t &mat_group,
        [[maybe_unused]] const index_t &index,
        [[maybe_unused]] const point2_t &loc,
        [[maybe_unused]] const scalar_t cos_inc_angle) const {

        using material_t = typename mat_group_t::value_type;

        constexpr auto inv{detail::invalid_value<scalar_t>()};

        // Access homogeneous surface material or material maps
        if constexpr ((detail::is_hom_material_v<material_t> &&
                       !std::is_same_v<material_t, material<scalar_t>>) ||
                      detail::is_material_map_v<material_t>) {

            // Slab or rod
            const auto mat =
                detail::material_accessor::get(mat_group, index, loc);

            // Empty material can occur in material maps, skip it
            if (!mat) {
                // Set the pathlength and thickness to zero so that they
                // are not counted
                return std::tuple(scalar_t{0}, scalar_t{0}, inv, inv);
            }

            const scalar_t seg{mat.path_segment(cos_inc_angle, loc[0])};
            const scalar_t t{mat.thickness()};
            const scalar_t mat_X0{mat.get_material().X0()};
            const scalar_t mat_L0{mat.get_material().L0()};

            return std::tuple(seg, t, mat_X0, mat_L0);
        } else {
            return std::tuple(inv, inv, inv, inv);
        }
    }
};

/// @brief Actor that collects all material encountered by a track during
///        navigation
///
/// The material is scaled with either the slab thickness or pathlength through
/// the material.
template <typename scalar_t>
struct material_tracer : detray::actor {

    struct state {
        material_record<scalar_t> mat_record{};
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state &tracer, const propagator_state_t &prop_state) const {

        const auto &navigation = prop_state._navigation;

        // Only count material if navigator encountered it
        if (!navigation.encountered_sf_material()) {
            return;
        }

        // For now use default context
        typename propagator_state_t::detector_type::geometry_context gctx{};

        // Current surface
        const auto sf = navigation.get_surface();

        // Bound position of track on current surface
        auto &bound_params = prop_state._stepping._bound_params;
        const auto loc_pos = bound_params.bound_local();
        const auto glob_dir = bound_params.dir();

        // Fetch the material parameters and pathlength through the material
        const auto [seg, t, mx0, ml0] =
            sf.template visit_material<get_material_params>(
                loc_pos, sf.cos_angle(gctx, glob_dir, loc_pos));

        // Fill the material record

        // Record the initial track direction
        if (detray::detail::is_invalid_value(tracer.mat_record.eta) &&
            detray::detail::is_invalid_value(tracer.mat_record.phi)) {
            tracer.mat_record.eta = getter::eta(glob_dir);
            tracer.mat_record.phi = getter::phi(glob_dir);
        }

        if (mx0 > 0.f) {
            tracer.mat_record.sX0 += seg / mx0;
            tracer.mat_record.tX0 += t / mx0;
        }
        if (ml0 > 0.f) {
            tracer.mat_record.sL0 += seg / ml0;
            tracer.mat_record.tL0 += t / ml0;
        }
    }
};

/// Write the accumulated material of a track from @param mat_records to a csv
/// file to the path @param mat_file_name
template <typename scalar_t>
auto write_material(const std::string &mat_file_name,
                    const dvector<material_record<scalar_t>> &mat_records) {

    const auto file_path = std::filesystem::path{mat_file_name};
    assert(file_path.extension() == ".csv");

    // Make sure path to file exists
    io::create_path(file_path.parent_path());

    detray::io::file_handle outfile{
        mat_file_name, std::ios::out | std::ios::binary | std::ios::trunc};
    *outfile << "eta,phi,mat_sX0,mat_sL0,mat_tX0,mat_tL0" << std::endl;

    for (const auto &rec : mat_records) {
        *outfile << rec.eta << "," << rec.phi << "," << rec.sX0 << ","
                 << rec.sL0 << "," << rec.tX0 << "," << rec.tL0 << std::endl;
    }
}

}  // namespace detray::material_validator
