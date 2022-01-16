/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/core/detector.hpp"
#include "tests/common/tools/detector_registry.hpp"

namespace detray {

template <template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple,
          template <typename...> class vector_type = dvector,
          template <typename...> class jagged_vector_type = djagged_vector>
auto create_toy_geometry(vecmem::memory_resource& resource) {

    // detector type
    using detector_t = detector<detector_registry::toy_detector, array_type,
                                tuple_type, vector_type, jagged_vector_type>;

    // sub-geometry components type
    using edge_links = typename detector_t::edge_type;
    using surface = typename detector_t::surface_type;
    using surfaces_regular_circular_grid =
        typename detector_t::surfaces_regular_circular_grid;
    using mask_container = typename detector_t::mask_container;
    using transform_store = typename detector_t::transform_store;
    using transform_container = typename detector_t::transform_container;
    using surface_container = typename detector_t::surface_filling_container;

    /** Function that adds a disc portal.
     */
    auto add_disc_portal =
        [](const dindex volume_id, typename transform_store::context& ctx,
           surface_container& surfaces, mask_container& masks,
           transform_container& transforms, const scalar min_r,
           const scalar max_r, const scalar half_z, const edge_links edge) {
            // translation
            point3 tsl{0., 0., half_z};

            // add transform
            transforms[detector_t::e_portal_ring2].emplace_back(ctx, tsl);

            // add mask
            masks.template add_mask<detector_t::e_portal_ring2>(min_r, max_r);
            masks.template group<detector_t::e_portal_ring2>().back().links() =
                edge;

            // create surface
            surface surf(
                transforms[detector_t::e_portal_ring2].size(ctx) - 1,
                {detector_t::e_portal_ring2,
                 masks.template group<detector_t::e_portal_ring2>().size() - 1},
                volume_id, dindex_invalid);

            surf.set_edge(edge);

            // add surface
            surfaces[detector_t::e_portal_ring2].push_back(surf);
        };

    /** Function that adds a disc portal.
     */
    auto add_cylinder_portal = [](const dindex volume_id,
                                  typename transform_store::context& ctx,
                                  surface_container& surfaces,
                                  mask_container& masks,
                                  transform_container& transforms,
                                  const scalar r, const scalar half_z,
                                  const edge_links edge) {
        // translation
        point3 tsl{0., 0., 0};

        // add transform
        transforms[detector_t::e_portal_cylinder3].emplace_back(ctx, tsl);

        // add mask
        masks.template add_mask<detector_t::e_portal_cylinder3>(r, -half_z,
                                                                half_z);
        masks.template group<detector_t::e_portal_cylinder3>().back().links() =
            edge;

        // create surface
        surface surf(
            transforms[detector_t::e_portal_cylinder3].size(ctx) - 1,
            {detector_t::e_portal_cylinder3,
             masks.template group<detector_t::e_portal_cylinder3>().size() - 1},
            volume_id, dindex_invalid);

        surf.set_edge(edge);

        // add surface
        surfaces[detector_t::e_portal_cylinder3].push_back(surf);
    };

    /** Function that creates the modules
     */
    auto create_modules =
        [&](const dindex volume_id, const dindex invalid_value,
            typename transform_store::context& ctx, surface_container& surfaces,
            mask_container& masks, transform_container& transforms,
            typename detector_t::surfaces_regular_circular_grid& surfaces_grid,
            const scalar m_half_x = 8.4, const scalar m_half_y = 36.,
            const scalar m_tilt_phi = 0.145, const scalar layer_r = 32.,
            const scalar radial_stagger = 2., const scalar l_overlap = 5.,
            const std::pair<int, int> binning = {16, 14}) {
            // surface grid bins
            int n_phi_bins = binning.first;
            int n_z_bins = binning.second;

            // module positions
            detray::dvector<point3> m_centers;
            m_centers.reserve(n_phi_bins * n_z_bins);

            // prep work
            scalar pi{static_cast<scalar>(M_PI)};
            scalar phi_step = scalar{2} * pi / (n_phi_bins);
            scalar min_phi = -pi + scalar{0.5} * phi_step;
            scalar max_phi = min_phi + n_phi_bins * phi_step;

            scalar z_start = scalar{-0.5} * (n_z_bins - 1) *
                             (scalar{2} * m_half_y - l_overlap);
            scalar z_step = scalar{2} * std::abs(z_start) / (n_z_bins - 1);
            scalar z_end = z_start + n_z_bins * z_step;

            // add surface grid
            typename detector_t::surfaces_circular_axis phi_axis(
                n_phi_bins, min_phi - phi_step * 0.5, max_phi - phi_step * 0.5,
                resource);
            typename detector_t::surfaces_regular_axis z_axis(
                n_z_bins, z_start - z_step * 0.5, z_end - z_step * 0.5,
                resource);

            surfaces_grid = typename detector_t::surfaces_regular_circular_grid(
                z_axis, phi_axis, resource);

            // loop over the bins
            for (size_t z_bin = 0; z_bin < size_t(n_z_bins); ++z_bin) {
                // prepare z and r
                scalar m_z = z_start + z_bin * z_step;
                scalar m_r = (z_bin % 2) != 0u
                                 ? layer_r - scalar{0.5} * radial_stagger
                                 : layer_r + scalar{0.5} * radial_stagger;
                for (size_t phiBin = 0; phiBin < size_t(n_phi_bins); ++phiBin) {
                    // calculate the current phi value
                    scalar m_phi = min_phi + phiBin * phi_step;
                    m_centers.push_back(point3{m_r * std::cos(m_phi),
                                               m_r * std::sin(m_phi), m_z});
                }
            }

            // Create surfaces
            for (auto& m_center : m_centers) {

                // Build the transform
                // The local phi
                scalar m_phi = algebra::getter::phi(m_center);
                // Local z axis is the normal vector
                vector3 m_local_z{std::cos(m_phi + m_tilt_phi),
                                  std::sin(m_phi + m_tilt_phi), 0.};
                // Local x axis the normal to local y,z
                vector3 m_local_x{-std::sin(m_phi + m_tilt_phi),
                                  std::cos(m_phi + m_tilt_phi), 0.};

                // add transform
                transforms[detector_t::e_rectangle2].emplace_back(
                    ctx, m_center, m_local_z, m_local_x);

                // add mask
                masks.template add_mask<detector_t::e_rectangle2>(m_half_x,
                                                                  m_half_y);
                masks.template group<detector_t::e_rectangle2>()
                    .back()
                    .links() = {volume_id, dindex_invalid};

                // create surface
                surface surf(
                    transforms[detector_t::e_rectangle2].size(ctx) - 1,
                    {detector_t::e_rectangle2,
                     masks.template group<detector_t::e_rectangle2>().size() -
                         1},
                    volume_id, dindex_invalid);

                surf.set_edge({volume_id, invalid_value});

                surf.set_grid_status(true);
                // add surface to surface container
                surfaces[detector_t::e_rectangle2].push_back(surf);
            }
        };

    /*********************************
     *  Let's build detector \(^0^)/ *
     *  Let's build detector \(^0^)/ *
     *  Let's build detector \(^0^)/ *
     *********************************/

    /**
     * Create detector and set some detector parameters
     */

    // create detector
    detector_t det(resource);

    // context objects
    typename transform_store::context ctx0;

    // detector parameters
    const scalar detector_half_z = 500.;
    const scalar beampipe_r = 27.;
    const scalar first_layer_outer_r = 38.;
    const scalar second_layer_inner_r = 64.;
    const scalar second_layer_outer_r = 80.;
    const dindex leaving_world = dindex_invalid;
    const dindex inv_sf_finder = dindex_invalid;

    /**
     * Create volumes
     */

    // beampipe volue -- volume ID = 0
    det.new_volume(
        {0, beampipe_r, -1 * detector_half_z, detector_half_z, -M_PI, M_PI});

    // first layer volume -- volume ID = 1
    det.new_volume({beampipe_r, first_layer_outer_r, -1 * detector_half_z,
                    detector_half_z, -M_PI, M_PI});

    // gap volume -- volume ID = 2
    det.new_volume({first_layer_outer_r, second_layer_inner_r,
                    -1 * detector_half_z, detector_half_z, -M_PI, M_PI});

    // second layer volume -- volume ID = 3
    det.new_volume({second_layer_inner_r, second_layer_outer_r,
                    -1 * detector_half_z, detector_half_z, -M_PI, M_PI});

    auto& vol0 = det.volume_by_index(0);
    auto& vol1 = det.volume_by_index(1);
    auto& vol2 = det.volume_by_index(2);
    auto& vol3 = det.volume_by_index(3);

    /**
     * Fill beampipe volume -- volume ID = 0
     */

    surface_container vol0_surfaces = {};
    mask_container vol0_masks(resource);
    transform_container vol0_transforms = {};

    // disc portal at negative side
    add_disc_portal(vol0.index(), ctx0, vol0_surfaces, vol0_masks,
                    vol0_transforms, 0, beampipe_r, -detector_half_z,
                    {leaving_world, inv_sf_finder});

    // disc_portal at positive side
    add_disc_portal(vol0.index(), ctx0, vol0_surfaces, vol0_masks,
                    vol0_transforms, 0, beampipe_r, detector_half_z,
                    {leaving_world, inv_sf_finder});

    // beampipe cylinder potal -> linked to first layer (vol1)
    add_cylinder_portal(vol0.index(), ctx0, vol0_surfaces, vol0_masks,
                        vol0_transforms, beampipe_r, detector_half_z,
                        {vol1.index(), inv_sf_finder});

    // Add all objects to detector
    det.add_objects(ctx0, vol0, vol0_surfaces, vol0_masks, vol0_transforms);

    /**
     * Fill the first layer volume -- volume ID = 1
     */

    surface_container vol1_surfaces = {};
    mask_container vol1_masks(resource);
    transform_container vol1_transforms = {};
    surfaces_regular_circular_grid vol1_surfaces_grid(resource);

    // disc portal at negative side
    add_disc_portal(vol1.index(), ctx0, vol1_surfaces, vol1_masks,
                    vol1_transforms, beampipe_r, first_layer_outer_r,
                    -detector_half_z, {leaving_world, inv_sf_finder});

    // disc portal at positive side
    add_disc_portal(vol1.index(), ctx0, vol1_surfaces, vol1_masks,
                    vol1_transforms, beampipe_r, first_layer_outer_r,
                    detector_half_z, {leaving_world, inv_sf_finder});

    // first layer inner cylinder potal -> linked to beampipe (vol0)
    add_cylinder_portal(vol1.index(), ctx0, vol1_surfaces, vol1_masks,
                        vol1_transforms, beampipe_r, detector_half_z,
                        {vol0.index(), inv_sf_finder});

    // first layer outer cylinder potal -> linked to gap volume (vol2)
    add_cylinder_portal(vol1.index(), ctx0, vol1_surfaces, vol1_masks,
                        vol1_transforms, first_layer_outer_r, detector_half_z,
                        {vol2.index(), inv_sf_finder});

    // create modules for the first layer
    create_modules(vol1.index(), inv_sf_finder, ctx0, vol1_surfaces, vol1_masks,
                   vol1_transforms, vol1_surfaces_grid, 8.4, 36., 0.145, 32.,
                   2., 5., {16, 14});

    // Add all objects to detector
    det.add_objects(ctx0, vol1, vol1_surfaces, vol1_masks, vol1_transforms);
    det.add_surfaces_grid(ctx0, vol1, vol1_surfaces_grid);

    /**
     * Fill the gap volume -- volume ID = 2
     */

    surface_container vol2_surfaces = {};
    mask_container vol2_masks(resource);
    transform_container vol2_transforms = {};

    // disc portal at negative side
    add_disc_portal(vol2.index(), ctx0, vol2_surfaces, vol2_masks,
                    vol2_transforms, first_layer_outer_r, second_layer_inner_r,
                    -detector_half_z, {leaving_world, inv_sf_finder});

    // disc portal at positive side
    add_disc_portal(vol2.index(), ctx0, vol2_surfaces, vol2_masks,
                    vol2_transforms, first_layer_outer_r, second_layer_inner_r,
                    detector_half_z, {leaving_world, inv_sf_finder});

    // gap volume inner cylinderical portal -> linked to first layer (vol1)
    add_cylinder_portal(vol2.index(), ctx0, vol2_surfaces, vol2_masks,
                        vol2_transforms, first_layer_outer_r, detector_half_z,
                        {vol1.index(), inv_sf_finder});

    // gap volume outer cylinderical portal -> linked to second layer (vol2)
    add_cylinder_portal(vol2.index(), ctx0, vol2_surfaces, vol2_masks,
                        vol2_transforms, second_layer_inner_r, detector_half_z,
                        {vol3.index(), inv_sf_finder});

    // Add all objects to detector
    det.add_objects(ctx0, vol2, vol2_surfaces, vol2_masks, vol2_transforms);

    /**
     * Fill the second -- volume ID = 3
     */

    surface_container vol3_surfaces = {};
    mask_container vol3_masks(resource);
    transform_container vol3_transforms = {};
    surfaces_regular_circular_grid vol3_surfaces_grid(resource);

    // disc portal at negative side
    add_disc_portal(vol3.index(), ctx0, vol3_surfaces, vol3_masks,
                    vol3_transforms, second_layer_inner_r, second_layer_outer_r,
                    -detector_half_z, {leaving_world, inv_sf_finder});

    // disc portal at positive side
    add_disc_portal(vol3.index(), ctx0, vol3_surfaces, vol3_masks,
                    vol3_transforms, second_layer_inner_r, second_layer_outer_r,
                    detector_half_z, {leaving_world, inv_sf_finder});

    // second layer inner cylinderical portal -> linked to gap volume (vol2)
    add_cylinder_portal(vol3.index(), ctx0, vol3_surfaces, vol3_masks,
                        vol3_transforms, second_layer_inner_r, detector_half_z,
                        {vol2.index(), inv_sf_finder});

    // second layer outer cylinderical portal -> linked to the end of world
    add_cylinder_portal(vol3.index(), ctx0, vol3_surfaces, vol3_masks,
                        vol3_transforms, second_layer_outer_r, detector_half_z,
                        {leaving_world, inv_sf_finder});

    // create modules for the second layer
    create_modules(vol3.index(), inv_sf_finder, ctx0, vol3_surfaces, vol3_masks,
                   vol3_transforms, vol3_surfaces_grid, 8.4, 36., 0.145, 72.,
                   2., 5., {32, 14});

    // Add all objects to detector
    det.add_objects(ctx0, vol3, vol3_surfaces, vol3_masks, vol3_transforms);
    det.add_surfaces_grid(ctx0, vol3, vol3_surfaces_grid);

    // return the detector
    return det;
}

}  // namespace detray