/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/rectangle2D.hpp"
#include "detray/tools/surface_factory_interface.hpp"

// System include(s)
#include <cassert>
#include <limits>

namespace detray {

/// @brief Generates a number of surfaces along a given direction
///
/// @tparam detector_t the type of detector the volume belongs to.
template <typename detector_t, typename mask_shape_t = rectangle2D<>,
          typename trajectory_t = detail::ray<typename detector_t::transform3>>
class telescope_generator final : public surface_factory_interface<detector_t> {

    using scalar_t = typename detector_t::scalar_type;
    using transform3_t = typename detector_t::transform3;
    using point3_t = typename detector_t::point3;
    using vector3_t = typename detector_t::vector3;

    public:
    /// Build a surface at with extent given in @param boundaries at every
    /// position in @param positions along the pilot-track @param traj.
    DETRAY_HOST
    telescope_generator(
        std::vector<scalar_t> positions,
        std::array<scalar_t, mask_shape_t::boundaries::e_size> boundaries,
        trajectory_t traj)
        : m_traj{traj}, m_positions{positions}, m_boundaries{boundaries} {}

    /// Infer the sensitive surface placement from the telescope @param length
    /// if no concrete positions were given.
    /// @param n_surfaces number of surfaces to be generated
    /// @param boundaries mask boundaries of the surfaces
    /// @param traj pilot track along which to build the telescope
    DETRAY_HOST
    telescope_generator(
        scalar_t length, std::size_t n_surfaces,
        std::array<scalar_t, mask_shape_t::boundaries::e_size> boundaries,
        trajectory_t traj)
        : m_traj{traj}, m_positions{}, m_boundaries{boundaries} {
        scalar_t pos{0.f};
        scalar_t dist{n_surfaces > 1u
                          ? length / static_cast<scalar_t>(n_surfaces - 1u)
                          : 0.f};
        for (std::size_t i = 0u; i < n_surfaces; ++i) {
            m_positions.push_back(pos);
            pos += dist;
        }
    }

    /// @returns the number of surfaces this factory will produce
    DETRAY_HOST
    auto size() const -> dindex override {
        return static_cast<dindex>(m_positions.size());
    }

    /// Clear the positions and boundaries of the surfaces.
    DETRAY_HOST
    void clear() override {
        m_positions.clear();
        m_boundaries = {};
    };

    /// This is a surface generator, no external surface data needed
    /// @{
    DETRAY_HOST
    void push_back(surface_data<detector_t> &&) override { /*Do nothing*/
    }
    DETRAY_HOST
    auto push_back(std::vector<surface_data<detector_t>> &&)
        -> void override { /*Do nothing*/
    }
    /// @}

    /// Create a surface telescope.
    ///
    /// @param volume the volume the portals need to be added to.
    /// @param surfaces the surface collection to wrap and to add the portals to
    /// @param transforms the transforms of the surfaces.
    /// @param masks the masks of the surfaces.
    /// @param ctx the geometry context (not needed for portals).
    DETRAY_HOST
    auto operator()(typename detector_t::volume_type &volume,
                    typename detector_t::surface_container &surfaces,
                    typename detector_t::transform_container &transforms,
                    typename detector_t::mask_container &masks,
                    typename detector_t::geometry_context ctx = {}) const
        -> dindex_range override {

        using surface_t = typename detector_t::surface_type;
        using nav_link_t = typename surface_t::navigation_link;
        using mask_link_t = typename surface_t::mask_link;
        using material_link_t = typename surface_t::material_link;

        const dindex surfaces_offset{static_cast<dindex>(surfaces.size())};

        // The type id of the surface mask shape
        constexpr auto mask_id{detector_t::mask_container::template get_id<
            mask<mask_shape_t>>::value};

        // The material will be added in a later step
        constexpr auto no_material{surface_t::material_id::e_none};

        auto volume_idx{volume.index()};

        // Create the module centers
        const auto mod_placements = module_positions(m_traj, m_positions);

        // Create geometry data
        for (const auto &mod_placement : mod_placements) {

            auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

            // Surfaces with the linking into the local containers
            mask_link_t mask_link{mask_id, masks.template size<mask_id>()};
            material_link_t material_link{
                no_material,
                detail::invalid_value<typename material_link_t::index_type>()};

            const auto trf_index = transforms.size(ctx);
            surfaces.emplace_back(trf_index, mask_link, material_link,
                                  volume_idx, dindex_invalid,
                                  surface_id::e_sensitive);

            // The rectangle bounds for this module
            masks.template emplace_back<mask_id>(empty_context{}, m_boundaries,
                                                 mask_volume_link);

            // Build the transform
            // Local z axis is the global normal vector
            vector3_t m_local_z = algebra::vector::normalize(mod_placement.dir);

            // Project onto the weakest direction component of the normal vector
            vector3_t e_i{0.f, 0.f, 0.f};
            auto min{std::numeric_limits<scalar_t>::max()};
            auto i{std::numeric_limits<unsigned int>::max()};
            for (unsigned int k = 0u; k < 3u; ++k) {
                if (m_local_z[k] < min) {
                    min = m_local_z[k];
                    i = k;
                }
            }
            e_i[i] = 1.f;
            vector3_t proj = algebra::vector::dot(m_local_z, e_i) * m_local_z;
            // Local x axis is the normal to local y,z
            vector3_t m_local_x = algebra::vector::normalize(e_i - proj);

            // Create the global-to-local transform of the module
            transforms.emplace_back(ctx, mod_placement.pos, m_local_z,
                                    m_local_x);
        }

        return {surfaces_offset, static_cast<dindex>(surfaces.size())};
    }

    private:
    /// Where and how to place the telescope modules.
    struct module_placement {
        /// Module position
        point3_t pos;
        /// Module normal
        vector3_t dir;
    };

    /// Helper method for positioning the surfaces.
    ///
    /// @param traj pilot trajectory along which the modules should be placed.
    /// @param steps lengths along the trajectory where surfaces should be
    ///              placed.
    ///
    /// @return a vector of the @c module_placements along the trajectory.
    inline std::vector<module_placement> module_positions(
        const trajectory_t &traj, const std::vector<scalar_t> &steps) const {

        // create and fill the module placements
        std::vector<module_placement> placements;
        placements.reserve(steps.size());

        for (const auto s : steps) {
            placements.push_back({traj.pos(s), traj.dir(s)});
        }

        return placements;
    }

    /// "pilot-track" along which to place the surfaces
    trajectory_t m_traj;
    /// Positions of the surfaces in the telescope along the pilot track
    std::vector<scalar_t> m_positions;
    /// The boundary values for the surface mask
    std::array<scalar_t, mask_shape_t::boundaries::e_size> m_boundaries;
};

}  // namespace detray
