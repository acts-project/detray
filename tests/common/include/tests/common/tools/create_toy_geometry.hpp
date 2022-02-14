/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <climits>
#include <stdexcept>
#include <type_traits>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/detector.hpp"
#include "tests/common/tools/detector_metadata.hpp"

namespace detray {

namespace {

/** Function that adds a cylinder portal.
 *
 * @tparam cylinder_id default cylinder id
 *
 * @param volume_id volume the portal should be added to
 * @param ctx geometric context
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param r raius of the cylinder
 * @param lower_z lower extend in z
 * @param upper_z upper extend in z
 * @param edge link to next volume and surfaces finder
 */
template <typename context_t, typename surface_container_t,
          typename mask_container_t, typename transform_container_t,
          typename edge_links>
inline void add_cylinder_surface(const dindex volume_id, context_t &ctx,
                                 surface_container_t &surfaces,
                                 mask_container_t &masks,
                                 transform_container_t &transforms,
                                 const scalar r, const scalar lower_z,
                                 const scalar upper_z, const edge_links edge) {
    using surface_t = typename surface_container_t::value_type::value_type;
    using mask_defs = typename surface_t::mask_defs;
    using mask_link_t = typename mask_container_t::link_type;

    constexpr auto cylinder_id = mask_defs::id::e_portal_cylinder3;

    const scalar min_z = std::min(lower_z, upper_z);
    const scalar max_z = std::max(lower_z, upper_z);

    // translation
    point3 tsl{0., 0., 0};

    // add transform and masks
    transforms[cylinder_id].emplace_back(ctx, tsl);
    masks.template add_mask<cylinder_id>(r, min_z, max_z, edge);

    // add surface
    mask_link_t mask_link{cylinder_id, masks.template size<cylinder_id>() - 1};
    const bool is_portal = std::get<0>(edge) != volume_id;
    surfaces[cylinder_id].emplace_back(transforms[cylinder_id].size(ctx) - 1,
                                       mask_link, volume_id, dindex_invalid,
                                       is_portal);
}

/** Function that adds a disc portal.
 *
 * @tparam disc_id default disc id
 *
 * @param volume_id volume the portal should be added to
 * @param ctx geometric context
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param min_r lower radius of the disc
 * @param max_r upper radius of the disc
 * @param z z position of the disc
 * @param edge link to next volume and surfaces finder
 */
template <typename context_t, typename surface_container_t,
          typename mask_container_t, typename transform_container_t,
          typename edge_links>
inline void add_disc_surface(const dindex volume_id, context_t &ctx,
                             surface_container_t &surfaces,
                             mask_container_t &masks,
                             transform_container_t &transforms,
                             const scalar inner_r, const scalar outer_r,
                             const scalar z, const edge_links edge) {
    using surface_t = typename surface_container_t::value_type::value_type;
    using mask_defs = typename surface_t::mask_defs;
    using mask_link_t = typename mask_container_t::link_type;

    constexpr auto disc_id = mask_defs::id::e_portal_ring2;

    const scalar min_r = std::min(inner_r, outer_r);
    const scalar max_r = std::max(inner_r, outer_r);

    // translation
    point3 tsl{0., 0., z};

    // add transform and mask
    transforms[disc_id].emplace_back(ctx, tsl);
    masks.template add_mask<disc_id>(min_r, max_r, edge);

    // add surface
    mask_link_t mask_link{disc_id, masks.template size<disc_id>() - 1};
    const bool is_portal = std::get<0>(edge) != volume_id;
    surfaces[disc_id].emplace_back(transforms[disc_id].size(ctx) - 1, mask_link,
                                   volume_id, dindex_invalid, is_portal);
}

/** Function that adds a generic cylinder volume, using a factory for contained
 *  module surfaces.
 *
 * @tparam detector_t the detector type
 * @tparam factory_t type of module factory. Must be callable on containers.
 *
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param lay_inner_r inner radius of volume
 * @param lay_outer_r outer radius of volume
 * @param lay_neg_r lower extend of volume
 * @param lay_pos_r upper extend of volume
 * @param edges volume and surfaces finder links for the portals of the volume
 * @param module_factory functor that adds module surfaces to volume
 */

template <
    typename detector_t, typename factory_t,
    std::enable_if_t<
        std::is_invocable_v<factory_t, typename detector_t::context &,
                            typename detector_t::volume_type &,
                            typename detector_t::surface_filling_container &,
                            typename detector_t::mask_container &,
                            typename detector_t::transform_filling_container &>,
        bool> = true>
void create_cyl_volume(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::context &ctx, const scalar lay_inner_r,
    const scalar lay_outer_r, const scalar lay_neg_z, const scalar lay_pos_z,
    const std::vector<typename detector_t::surface_type::edge_type> &edges,
    factory_t &module_factory) {
    // volume bounds
    const scalar inner_r = std::min(lay_inner_r, lay_outer_r);
    const scalar outer_r = std::max(lay_inner_r, lay_outer_r);
    const scalar lower_z = std::min(lay_neg_z, lay_pos_z);
    const scalar upper_z = std::max(lay_neg_z, lay_pos_z);

    auto &cyl_volume =
        det.new_volume({inner_r, outer_r, lower_z, upper_z, -M_PI, M_PI});
    typename detector_t::surfaces_finder_type::surfaces_regular_circular_grid
        cyl_surfaces_grid(resource);

    // Add module surfaces to volume
    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_filling_container transforms = {resource};

    // fill the surfaces
    module_factory(ctx, cyl_volume, surfaces, masks, transforms);
    // create the surface grid
    module_factory(cyl_surfaces_grid, *det.resource());

    // negative and positive, inner and outer portal surface
    add_cylinder_surface(cyl_volume.index(), ctx, surfaces, masks, transforms,
                         inner_r, lower_z, upper_z, edges[0]);
    add_cylinder_surface(cyl_volume.index(), ctx, surfaces, masks, transforms,
                         outer_r, lower_z, upper_z, edges[1]);
    add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, transforms,
                     inner_r, outer_r, lower_z, edges[2]);
    add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, transforms,
                     inner_r, outer_r, upper_z, edges[3]);

    det.add_objects(ctx, cyl_volume, surfaces, masks, transforms);
    if (cyl_volume.get_grid_type() !=
        detector_t::volume_type::grid_type::e_no_grid) {
        det.add_surfaces_grid(ctx, cyl_volume, cyl_surfaces_grid);
    }
}

/** Helper function that creates a layer of rectangular barrel modules.
 *
 * @tparam rctangle_id default rectangle id
 *
 * @param ctx geometric context
 * @param volume_id volume the portal should be added to
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param cfg config struct for module creation
 */
template <typename context_t, typename volume_type,
          typename surface_container_t, typename mask_container_t,
          typename transform_container_t, typename config_t>
inline void create_barrel_modules(context_t &ctx, volume_type &vol,
                                  surface_container_t &surfaces,
                                  mask_container_t &masks,
                                  transform_container_t &transforms,
                                  config_t cfg) {
    using surface_t = typename surface_container_t::value_type::value_type;
    using mask_defs = typename surface_t::mask_defs;
    using edge_t = typename surface_t::edge_type;
    using mask_link_t = typename mask_container_t::link_type;

    constexpr auto rectangle_id = mask_defs::id::e_rectangle2;
    auto volume_id = vol.index();
    edge_t mask_edge{volume_id, dindex_invalid};

    vol.set_grid_type(volume_type::grid_type::e_z_phi_grid);

    // Create the module centers

    // surface grid bins
    int n_phi_bins = cfg.m_binning.first;
    int n_z_bins = cfg.m_binning.second;
    // module positions
    detray::dvector<point3> m_centers;
    m_centers.reserve(n_phi_bins * n_z_bins);

    // prep work
    // scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step = 2 * M_PI / (n_phi_bins);
    scalar min_phi = -M_PI + 0.5 * phi_step;

    auto z_axis_info = cfg.get_z_axis_info();
    auto &z_start = std::get<0>(z_axis_info);
    // auto &z_end = std::get<1>(z_axis_info);
    auto &z_step = std::get<2>(z_axis_info);

    // loop over the z bins
    for (size_t z_bin = 0; z_bin < size_t(n_z_bins); ++z_bin) {
        // prepare z and r
        scalar m_z = z_start + z_bin * z_step;
        scalar m_r = (z_bin % 2) != 0u
                         ? cfg.layer_r - 0.5 * cfg.m_radial_stagger
                         : cfg.layer_r + 0.5 * cfg.m_radial_stagger;
        for (size_t phiBin = 0; phiBin < size_t(n_phi_bins); ++phiBin) {
            // calculate the current phi value
            scalar m_phi = min_phi + phiBin * phi_step;
            m_centers.push_back(
                point3{m_r * std::cos(m_phi), m_r * std::sin(m_phi), m_z});
        }
    }

    // Create geometry data
    for (auto &m_center : m_centers) {

        // Surfaces with the linking into the local containers
        mask_link_t m_id = {rectangle_id, masks.template size<rectangle_id>()};
        const auto trf_index = transforms[rectangle_id].size(ctx);
        surfaces[rectangle_id].emplace_back(trf_index, m_id, volume_id,
                                            dindex_invalid, false);
        surfaces[rectangle_id].back().set_grid_status(true);

        // The rectangle bounds for this module
        masks.template add_mask<rectangle_id>(cfg.m_half_x, cfg.m_half_y,
                                              mask_edge);

        // Build the transform
        // The local phi
        scalar m_phi = algebra::getter::phi(m_center);
        // Local z axis is the normal vector
        vector3 m_local_z{std::cos(m_phi + cfg.m_tilt_phi),
                          std::sin(m_phi + cfg.m_tilt_phi), 0.};
        // Local x axis the normal to local y,z
        vector3 m_local_x{-std::sin(m_phi + cfg.m_tilt_phi),
                          std::cos(m_phi + cfg.m_tilt_phi), 0.};

        // Create the module transform
        transforms[rectangle_id].emplace_back(ctx, m_center, m_local_z,
                                              m_local_x);
    }
}

/** Helper function that creates a surface grid of rectangular barrel modules.
 *
 * @param surfaces_grid the grid to be created with proper axes
 * @param resource vecmem memory resource
 * @param cfg config struct for module creation
 */
template <typename surfaces_grid_t, typename config_t>
inline void create_barrel_grid(surfaces_grid_t &surfaces_grid,
                               vecmem::memory_resource &resource,
                               config_t cfg) {
    auto z_axis_info = cfg.get_z_axis_info();
    auto &z_start = std::get<0>(z_axis_info);
    auto &z_end = std::get<1>(z_axis_info);
    auto &z_step = std::get<2>(z_axis_info);

    // add surface grid
    typename surfaces_grid_t::axis_p0_type z_axis(
        cfg.m_binning.second, z_start - z_step * 0.5, z_end + z_step * 0.5,
        resource);
    typename surfaces_grid_t::axis_p1_type phi_axis(cfg.m_binning.first, -M_PI,
                                                    M_PI, resource);

    surfaces_grid = surfaces_grid_t(z_axis, phi_axis, resource);
}

/** Helper method for positioning of modules in an endcap ring
 *
 * @param z is the z position of the ring
 * @param radius is the ring radius
 * @param phi_stagger is the radial staggering along phi
 * @param phi_sub_stagger is the overlap of the modules
 * @param n_phi_bins is the number of bins in phi
 *
 * @return a vector of the module positions in a ring
 */
inline auto module_positions_ring(scalar z, scalar radius, scalar phi_stagger,
                                  scalar phi_sub_stagger, int n_phi_bins) {
    // create and fill the positions
    std::vector<vector3> r_positions;
    r_positions.reserve(n_phi_bins);

    // prep work
    // scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step = scalar{2} * M_PI / (n_phi_bins);
    scalar min_phi = -M_PI + 0.5 * phi_step;

    for (size_t iphi = 0; iphi < size_t(n_phi_bins); ++iphi) {
        // if we have a phi sub stagger presents
        scalar rzs = 0.;
        // phi stagger affects 0 vs 1, 2 vs 3 ... etc
        // -> only works if it is a %4
        // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
        if (phi_sub_stagger != 0. && !(n_phi_bins % 4)) {
            // switch sides
            if (!(iphi % 4)) {
                rzs = phi_sub_stagger;
            } else if (!((iphi + 1) % 4)) {
                rzs = -phi_sub_stagger;
            }
        }
        // the module phi
        scalar phi = min_phi + iphi * phi_step;
        // main z position depending on phi bin
        scalar rz = iphi % 2 ? z - 0.5 * phi_stagger : z + 0.5 * phi_stagger;
        r_positions.push_back(
            vector3{radius * std::cos(phi), radius * std::sin(phi), rz + rzs});
    }
    return r_positions;
}

/** Helper function that creates a layer of trapezoidal endcap modules.
 *
 * @tparam trapezoid_id default trapezoid id
 *
 * @param ctx geometric context
 * @param volume_id volume the portal should be added to
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param cfg config struct for module creation
 */
template <typename context_t, typename volume_type,
          typename surface_container_t, typename mask_container_t,
          typename transform_container_t, typename config_t>
void create_endcap_modules(context_t &ctx, volume_type &vol,
                           surface_container_t &surfaces,
                           mask_container_t &masks,
                           transform_container_t &transforms, config_t cfg) {
    using surface_t = typename surface_container_t::value_type::value_type;
    using mask_defs = typename surface_t::mask_defs;
    using edge_t = typename surface_t::edge_type;
    using mask_link_t = typename mask_container_t::link_type;

    constexpr auto trapezoid_id = mask_defs::id::e_trapezoid2;
    auto volume_id = vol.index();
    edge_t mask_edge{volume_id, dindex_invalid};

    vol.set_grid_type(volume_type::grid_type::e_r_phi_grid);

    // calculate the radii of the rings
    std::vector<scalar> radii;
    // calculate the radial borders
    // std::vector<scalar> radial_boarders;
    // the radial span of the disc
    scalar delta_r = cfg.outer_r - cfg.inner_r;

    // Only one ring
    if (cfg.disc_binning.size() == 1) {
        radii.push_back(scalar{0.5} * (cfg.inner_r + cfg.outer_r));
        // radial_boarders = {inner_r, outer_r};
    } else {
        // sum up the total length of the modules along r
        scalar tot_length = 0;
        for (auto &m_hlength : cfg.m_half_y) {
            tot_length += 2 * m_hlength + 0.5;
        }
        // now calculate the overlap (equal pay)
        scalar r_overlap = (tot_length - delta_r) / (cfg.m_half_y.size() - 1);
        // and now fill the radii and gaps
        scalar prev_r = cfg.inner_r;
        scalar prev_hl = 0.;
        scalar prev_ol = 0.;
        // remember the radial boarders
        // radial_boarders.push_back(inner_r);
        for (auto &m_hlength : cfg.m_half_y) {
            // calculate the radius
            radii.push_back(prev_r + prev_hl - prev_ol + m_hlength);
            prev_r = radii.back();
            prev_ol = r_overlap;
            prev_hl = m_hlength;
            // and register the radial boarder
            // radial_boarders.push_back(prev_r + scalar{2} * prev_hl -
            // scalar{0.5} * prev_ol);
        }
    }

    // now build the modules in every ring
    for (size_t ir = 0; ir < radii.size(); ++ir) {
        // generate the z value
        // convention inner ring is closer to origin : makes sense
        scalar rz =
            radii.size() == 1
                ? cfg.edc_position
                : (ir % 2 ? cfg.edc_position + scalar{0.5} * cfg.ring_stagger
                          : cfg.edc_position - scalar{0.5} * cfg.ring_stagger);
        // fill the ring module positions
        scalar ps_stagger =
            cfg.m_phi_sub_stagger.size() ? cfg.m_phi_sub_stagger[ir] : 0.;

        std::vector<point3> r_postitions =
            module_positions_ring(rz, radii[ir], cfg.m_phi_stagger[ir],
                                  ps_stagger, cfg.disc_binning[ir]);

        // Build the geometrical objects
        for (const auto &m_position : r_postitions) {
            // trapezoid mask
            mask_link_t mask_link{trapezoid_id,
                                  masks.template size<trapezoid_id>()};
            masks.template add_mask<trapezoid_id>(cfg.m_half_x_min_y[ir],
                                                  cfg.m_half_x_max_y[ir],
                                                  cfg.m_half_y[ir], mask_edge);

            // Surfaces with the linking into the local containers
            surfaces[trapezoid_id].emplace_back(
                transforms[trapezoid_id].size(ctx), mask_link, volume_id,
                dindex_invalid, false);
            surfaces[trapezoid_id].back().set_grid_status(true);

            // the module transform from the position
            scalar m_phi = algebra::getter::phi(m_position);
            // the center position of the modules
            point3 m_center{static_cast<scalar>(cfg.side) * m_position};
            // the rotation matrix of the module
            vector3 m_local_y{std::cos(m_phi), std::sin(m_phi), 0.};
            // take different axis to have the same readout direction
            vector3 m_local_z{0., 0., cfg.side * scalar{1.}};
            vector3 m_local_x = algebra::vector::cross(m_local_y, m_local_z);

            // Create the module transform
            transforms[trapezoid_id].emplace_back(ctx, m_center, m_local_z,
                                                  m_local_x);
        }
    }
}

/** Helper function that creates a surface grid of trapezoidal endcap modules.
 *
 * @param surfaces_grid the grid to be created with proper axes
 * @param resource vecmem memory resource
 * @param cfg config struct for module creation
 */
template <typename surfaces_grid_t, typename config_t>
inline void create_endcap_grid(surfaces_grid_t &surfaces_grid,
                               vecmem::memory_resource &resource,
                               config_t cfg) {
    // add surface grid
    // TODO: WHat is the proper value of n_phi_bins?
    typename surfaces_grid_t::axis_p0_type r_axis(
        cfg.disc_binning.size(), cfg.inner_r, cfg.outer_r, resource);
    typename surfaces_grid_t::axis_p1_type phi_axis(cfg.disc_binning.front(),
                                                    -M_PI, M_PI, resource);

    surfaces_grid = surfaces_grid_t(r_axis, phi_axis, resource);
}

/** Helper method for creating a beampipe with enclosing volume.
 *
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_layers number of endcap layers that contain modules
 * @param edc_lay_sizes extend of the outer cylinder portal surfaces in z
 * @param beampipe_vol_size inner and outer radious of the beampipe volume
 * @param beampipe_r radius of the beampipe surface
 * @param brl_half_z half length of the barrel region in z
 */
template <typename detector_t>
inline void add_beampipe(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::context &ctx, const std::size_t n_edc_layers,
    const std::size_t n_brl_layers,
    const std::vector<std::pair<scalar, scalar>> &edc_lay_sizes,
    const std::pair<scalar, scalar> &beampipe_vol_size, const scalar beampipe_r,
    const scalar brl_half_z, const scalar edc_inner_r) {

    const dindex inv_sf_finder = dindex_invalid;
    const dindex leaving_world = dindex_invalid;

    scalar max_z =
        n_edc_layers <= 0 ? brl_half_z : edc_lay_sizes[n_edc_layers - 1].second;
    scalar min_z = -max_z;

    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_filling_container transforms = {resource};

    auto &beampipe =
        det.new_volume({beampipe_vol_size.first, beampipe_vol_size.second,
                        min_z, max_z, -M_PI, M_PI});
    const auto beampipe_idx = beampipe.index();

    // This is the beampipe surface
    typename detector_t::surface_type::edge_type edge = {beampipe_idx,
                                                         inv_sf_finder};
    add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, transforms,
                         beampipe_r, min_z, max_z, edge);

    // Get vol sizes in z, including for gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {brl_half_z, edc_lay_sizes[0].first}};
    for (std::size_t i = 0; i < n_edc_layers; ++i) {
        vol_sizes.emplace_back(edc_lay_sizes[i].first, edc_lay_sizes[i].second);
        vol_sizes.emplace_back(edc_lay_sizes[i].second,
                               edc_lay_sizes[i + 1].first);
    }
    vol_sizes.pop_back();

    // negative endcap portals
    unsigned int volume_link = beampipe_idx;
    for (int i = vol_sizes.size() - 1; i >= 0; --i) {
        edge = {++volume_link, inv_sf_finder};
        add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, transforms,
                             edc_inner_r, -vol_sizes[i].second,
                             -vol_sizes[i].first, edge);
    }
    // barrel portals
    if (n_brl_layers <= 0) {
        edge = {leaving_world, inv_sf_finder};
    } else {
        edge = {volume_link + 1, inv_sf_finder};
    }
    add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, transforms,
                         edc_inner_r, -brl_half_z, brl_half_z, edge);

    // positive endcap portals
    volume_link += 7;
    for (std::size_t i = 0; i < vol_sizes.size(); ++i) {
        edge = {++volume_link, inv_sf_finder};
        add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, transforms,
                             edc_inner_r, vol_sizes[i].second,
                             vol_sizes[i].first, edge);
    }

    // disc portals
    edge = {leaving_world, inv_sf_finder};
    add_disc_surface(beampipe_idx, ctx, surfaces, masks, transforms,
                     beampipe_vol_size.first, beampipe_vol_size.second, min_z,
                     edge);
    add_disc_surface(beampipe_idx, ctx, surfaces, masks, transforms,
                     beampipe_vol_size.first, beampipe_vol_size.second, max_z,
                     edge);

    det.add_objects(ctx, beampipe, surfaces, masks, transforms);
}

/** Helper method for creating a connecting gap volume between endcap and barrel
 *
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_brl_layers number of barrel layers that contain modules
 * @param beampipe_idx index of the beampipe volume
 * @param brl_lay_sizes extend of the disc portal surfaces in r
 * @param edc_inner_r inner radius of the gap volume
 * @param edc_outer_r outer radius of the gap volume
 * @param gap_neg_z lower extend of the gap volume in z
 * @param gap_pos_z upper extend of the gap volume in z
 * @param brl_vol_idx index of the first barrel volume (innermost layer)
 * @param edc_vol_idx index of the bordering endcap volume
 */
template <typename detector_t>
inline void add_endcap_barrel_connection(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::context &ctx, const int side,
    const unsigned int n_brl_layers, const dindex beampipe_idx,
    const std::vector<std::pair<scalar, scalar>> &brl_lay_sizes,
    const scalar edc_inner_r, const scalar edc_outer_r,
    const scalar gap_lower_z, const scalar gap_upper_z, dindex brl_vol_idx,
    const dindex edc_vol_idx) {
    const scalar min_z = std::min(side * gap_lower_z, side * gap_upper_z);
    const scalar max_z = std::max(side * gap_lower_z, side * gap_upper_z);
    scalar edc_disc_z = side < 0 ? min_z : max_z;
    scalar brl_disc_z = side < 0 ? max_z : min_z;

    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_filling_container transforms = {resource};

    auto &connector_gap =
        det.new_volume({edc_inner_r, edc_outer_r, min_z, max_z, -M_PI, M_PI});
    dindex connector_gap_idx = det.volumes().back().index();
    dindex leaving_world = dindex_invalid, inv_sf_finder = dindex_invalid;

    typename detector_t::surface_type::edge_type edge = {beampipe_idx,
                                                         inv_sf_finder};
    add_cylinder_surface(connector_gap_idx, ctx, surfaces, masks, transforms,
                         edc_inner_r, min_z, max_z, edge);
    edge = {leaving_world, inv_sf_finder};
    add_cylinder_surface(connector_gap_idx, ctx, surfaces, masks, transforms,
                         edc_outer_r, min_z, max_z, edge);
    edge = {edc_vol_idx, inv_sf_finder};
    add_disc_surface(connector_gap_idx, ctx, surfaces, masks, transforms,
                     edc_inner_r, edc_outer_r, edc_disc_z, edge);

    // Get vol sizes in z also for gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes;
    for (std::size_t i = 1; i <= n_brl_layers; ++i) {
        vol_sizes.emplace_back(brl_lay_sizes[i].first, brl_lay_sizes[i].second);
        vol_sizes.emplace_back(brl_lay_sizes[i].second,
                               brl_lay_sizes[i + 1].first);
    }

    edge = {brl_vol_idx, inv_sf_finder};
    for (std::size_t i = 0; i < 2 * n_brl_layers - 1; ++i) {
        edge = {brl_vol_idx++, inv_sf_finder};
        add_disc_surface(connector_gap_idx, ctx, surfaces, masks, transforms,
                         vol_sizes[i].first, vol_sizes[i].second, brl_disc_z,
                         edge);
    }

    det.add_objects(ctx, connector_gap, surfaces, masks, transforms);
}

/** Helper method for creating one of the two endcaps.
 *
 * @param det detector the subdetector should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_layers number of layers that contain modules
 * @param beampipe_idx index of the beampipe outermost volume
 * @param lay_sizes extend of the endcap layers in z direction
 * @param lay_positions position of the endcap layers in z direction
 * @param cfg config struct for module creation
 */
template <typename empty_vol_factory, typename edc_module_factory,
          typename detector_t, typename config_t>
void add_endcap_detector(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::context &ctx, std::size_t n_layers,
    dindex beampipe_idx,
    const std::vector<std::pair<scalar, scalar>> &lay_sizes,
    const std::vector<scalar> &lay_positions, config_t cfg) {

    // Generate consecutive linking between volumes (all edges for every vol.)
    using edge_t = typename detector_t::surface_type::edge_type;
    std::vector<std::vector<edge_t>> edges_vec;
    dindex leaving_world = dindex_invalid, inv_sf_finder = dindex_invalid;
    dindex first_vol_idx = det.volumes().back().index() + 1;
    dindex last_vol_idx = first_vol_idx + 2 * n_layers - 2;
    dindex prev_vol_idx = first_vol_idx - 1;
    dindex next_vol_idx = first_vol_idx + 1;

    for (int i = 0; i < 2 * static_cast<int>(n_layers) - 3; ++i) {
        edges_vec.push_back({{beampipe_idx, inv_sf_finder},
                             {leaving_world, inv_sf_finder},
                             {++prev_vol_idx, inv_sf_finder},
                             {++next_vol_idx, inv_sf_finder}});
    }
    // Edge of the world is flipped
    if (cfg.side < 0) {
        edges_vec.insert(edges_vec.begin(),
                         {{beampipe_idx, inv_sf_finder},
                          {leaving_world, inv_sf_finder},
                          {leaving_world, inv_sf_finder},
                          {first_vol_idx + 1, inv_sf_finder}});

        edges_vec.push_back({{beampipe_idx, inv_sf_finder},
                             {leaving_world, inv_sf_finder},
                             {last_vol_idx - 1, inv_sf_finder},
                             {last_vol_idx + 1, inv_sf_finder}});
    } else {
        // For n_layers=1 no extra gap layer is needed
        if (n_layers > 1) {
            edges_vec.insert(edges_vec.begin(),
                             {{beampipe_idx, inv_sf_finder},
                              {leaving_world, inv_sf_finder},
                              {first_vol_idx - 1, inv_sf_finder},
                              {first_vol_idx + 1, inv_sf_finder}});
        }

        edges_vec.push_back({{beampipe_idx, inv_sf_finder},
                             {leaving_world, inv_sf_finder},
                             {last_vol_idx - 1, inv_sf_finder},
                             {leaving_world, inv_sf_finder}});
    }

    // Get vol sizes in z, including gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {lay_sizes[0].first, lay_sizes[0].second}};
    for (std::size_t i = 1; i < n_layers; ++i) {
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i - 1].second);
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i].second);
    }

    edc_module_factory m_factory{cfg};
    empty_vol_factory empty_factory{};

    auto vol_size_itr = vol_sizes.begin();
    auto pos_itr = lay_positions.begin();
    // Reverse iteration for negative endcap
    if (cfg.side < 0) {
        std::advance(vol_size_itr, 2 * n_layers - 2);
        std::advance(pos_itr, n_layers - 1);
    }
    bool is_gap = true;
    for (std::size_t i = 0; i < 2 * n_layers - 1; ++i) {
        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {
            create_cyl_volume(det, resource, ctx, cfg.inner_r, cfg.outer_r,
                              cfg.side * (vol_size_itr + cfg.side * i)->first,
                              cfg.side * (vol_size_itr + cfg.side * i)->second,
                              edges_vec[i], empty_factory);

        } else {
            m_factory.cfg.edc_position = *(pos_itr + cfg.side * i / 2);
            create_cyl_volume(det, resource, ctx, cfg.inner_r, cfg.outer_r,
                              cfg.side * (vol_size_itr + cfg.side * i)->first,
                              cfg.side * (vol_size_itr + cfg.side * i)->second,
                              edges_vec[i], m_factory);
        }
    }
}

/** Helper method for creating the barrel section.
 *
 * @param det detector the subdetector should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_layers number of layers that contain modules
 * @param beampipe_idx index of the beampipe outermost volume
 * @param brl_half_z half length of the barrel section in z direction
 * @param lay_sizes extend of the barrel layers in r direction
 * @param lay_positions position of the barrel layers in r direction
 * @param cfg config struct for module creation
 */
template <typename empty_vol_factory, typename brl_module_factory,
          typename detector_t, typename config_t>
void add_barrel_detector(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::context &ctx, const unsigned int n_layers,
    dindex beampipe_idx, const scalar brl_half_z,
    const std::vector<std::pair<scalar, scalar>> &lay_sizes,
    const std::vector<scalar> &lay_positions,
    const std::vector<std::pair<int, int>> &m_binning, config_t cfg) {

    // Generate consecutive linking between volumes
    dindex leaving_world = dindex_invalid, inv_sf_finder = dindex_invalid;
    dindex first_vol_idx = det.volumes().back().index();
    dindex last_vol_idx = first_vol_idx + 2 * n_layers;
    dindex prev_vol_idx = first_vol_idx;
    dindex next_vol_idx = n_layers > 1 ? first_vol_idx + 2 : leaving_world;

    // Leave world, if no endcaps are present
    if (det.volumes().back().index() == 0) {
        first_vol_idx = leaving_world;
        last_vol_idx = leaving_world;
    }

    // First barrel layer is connected to the beampipe
    using edge_t = typename detector_t::surface_type::edge_type;
    std::vector<std::vector<edge_t>> edges_vec{{{beampipe_idx, inv_sf_finder},
                                                {next_vol_idx, inv_sf_finder},
                                                {first_vol_idx, inv_sf_finder},
                                                {last_vol_idx, inv_sf_finder}}};

    for (std::size_t i = 1; i < 2 * n_layers - 2; ++i) {
        edges_vec.push_back({{++prev_vol_idx, inv_sf_finder},
                             {++next_vol_idx, inv_sf_finder},
                             {first_vol_idx, inv_sf_finder},
                             {last_vol_idx, inv_sf_finder}});
    }
    // Last barrel layer leaves detector world
    edges_vec.push_back({{++prev_vol_idx, inv_sf_finder},
                         {leaving_world, inv_sf_finder},
                         {first_vol_idx, inv_sf_finder},
                         {last_vol_idx, inv_sf_finder}});

    // Get vol sizes in z, including gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {lay_sizes[1].first, lay_sizes[1].second}};
    for (std::size_t i = 2; i < n_layers + 1; ++i) {
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i - 1].second);
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i].second);
    }

    brl_module_factory m_factory{cfg};
    empty_vol_factory empty_factory{};
    bool is_gap = true;
    for (unsigned int i = 0; i < 2 * n_layers - 1; ++i) {
        unsigned int j = (i + 2) / 2;
        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {
            create_cyl_volume(det, resource, ctx, vol_sizes[i].first,
                              vol_sizes[i].second, -brl_half_z, brl_half_z,
                              edges_vec[i], empty_factory);
        } else {
            m_factory.cfg.m_binning = m_binning[j];
            m_factory.cfg.layer_r = lay_positions[j];
            create_cyl_volume(det, resource, ctx, vol_sizes[i].first,
                              vol_sizes[i].second, -brl_half_z, brl_half_z,
                              edges_vec[i], m_factory);
        }
    }
}

}  // namespace

/** Builds a detray geometry that contains the innermost tml layers. The number
 *  of barrel and endcap layers can be chosen, but all barrel layers should be
 *  present when an endcap detector is built to have the barrel region radius
 *  match the endcap diameter.
 *
 * @param n_brl_layers number of pixel barrel layer to build (max 4)
 * @param n_edc_layers number of pixel endcap discs to build (max 7)
 *
 * @returns a complete detector object
 */
template <template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple,
          template <typename...> class vector_t = dvector,
          template <typename...> class jagged_vector_t = djagged_vector>
auto create_toy_geometry(vecmem::memory_resource &resource,
                         std::size_t n_brl_layers = 4,
                         std::size_t n_edc_layers = 3) {

    // detector type
    using detector_t = detector<detector_registry::toy_detector, array_t,
                                tuple_t, vector_t, jagged_vector_t>;

    /** Leaving world */
    const dindex leaving_world = dindex_invalid;

    //
    // barrel
    //
    const scalar brl_half_z = 500.;
    const std::vector<scalar> brl_positions = {19., 32., 72., 116., 172.};
    const std::vector<std::pair<scalar, scalar>> brl_lay_sizes = {
        {0., 27.}, {27., 38.}, {64., 80.}, {108., 124.}, {164., 180.}};
    const std::vector<std::pair<int, int>> brl_binning = {
        {0., 0.}, {16, 14}, {32, 14}, {52, 14}, {78, 14}};
    // module parameters
    struct brl_m_config {
        scalar m_half_x = 8.4;
        scalar m_half_y = 36.;
        scalar m_tilt_phi = 0.14;  // 0.145;
        scalar layer_r = 32.;
        scalar m_radial_stagger = 0.5;  // 2.;
        scalar m_long_overlap = 2.;     // 5.;
        std::pair<int, int> m_binning = {16, 14};

        // return the first z position of module
        std::tuple<scalar, scalar, scalar> get_z_axis_info() {
            auto n_z_bins = m_binning.second;
            scalar z_start =
                -0.5 * (n_z_bins - 1) * (2 * m_half_y - m_long_overlap);
            scalar z_end = std::abs(z_start);
            scalar z_step = (z_end - z_start) / (n_z_bins - 1);

            return {z_start, z_end, z_step};
        }
    };

    //
    // endcaps
    //
    const std::vector<scalar> edc_positions = {600.,  700.,  820., 960.,
                                               1100., 1300., 1500.};
    const std::vector<std::pair<scalar, scalar>> edc_lay_sizes = {
        {595., 605.},   {695., 705.},   {815., 825.},  {955., 965.},
        {1095., 1105.}, {1295., 1305.}, {1495., 1505.}};
    // module params
    struct edc_m_config {
        int side = 1;
        scalar inner_r = 27.;
        scalar outer_r = 180.;
        scalar edc_position = 600.;
        scalar ring_stagger = 1.0;
        // Parameters for both rings of modules
        std::vector<scalar> m_phi_stagger = {4.0, 4.0};
        std::vector<scalar> m_phi_sub_stagger = {0.5, 0.};
        std::vector<size_t> disc_binning = {40, 68};
        std::vector<scalar> m_half_y = {36., 36.};
        std::vector<scalar> m_half_x_min_y = {8.4, 8.4};
        std::vector<scalar> m_half_x_max_y = {12.4, 12.4};
        std::vector<scalar> m_tilt = {0., 0.};
    };

    // Don't create modules in gap volume
    struct empty_vol_factory {
        void operator()(
            typename detector_t::context & /*ctx*/,
            typename detector_t::volume_type & /*volume*/,
            typename detector_t::surface_filling_container & /*surfaces*/,
            typename detector_t::mask_container & /*masks*/,
            typename detector_t::transform_filling_container & /*transforms*/) {
        }
        void operator()(typename detector_t::surfaces_finder_type::
                            surfaces_regular_circular_grid & /*surfaces_grid*/,
                        vecmem::memory_resource & /*resource*/) {}
    };

    // Fills volume with barrel layer
    struct brl_module_factory {
        brl_m_config cfg;

        void operator()(
            typename detector_t::context &ctx,
            typename detector_t::volume_type &volume,
            typename detector_t::surface_filling_container &surfaces,
            typename detector_t::mask_container &masks,
            typename detector_t::transform_filling_container &transforms) {
            create_barrel_modules(ctx, volume, surfaces, masks, transforms,
                                  cfg);
        }
        void operator()(typename detector_t::surfaces_finder_type::
                            surfaces_regular_circular_grid &surfaces_grid,
                        vecmem::memory_resource &resource) {
            create_barrel_grid(surfaces_grid, resource, cfg);
        }
    };

    // Fills volume with endcap rings
    struct edc_module_factory {
        edc_m_config cfg;

        void operator()(
            typename detector_t::context &ctx,
            typename detector_t::volume_type &volume,
            typename detector_t::surface_filling_container &surfaces,
            typename detector_t::mask_container &masks,
            typename detector_t::transform_filling_container &transforms) {
            create_endcap_modules(ctx, volume, surfaces, masks, transforms,
                                  cfg);
        }
        void operator()(typename detector_t::surfaces_finder_type::
                            surfaces_regular_circular_grid &surfaces_grid,
                        vecmem::memory_resource &resource) {
            create_endcap_grid(surfaces_grid, resource, cfg);
        }
    };

    // create empty detector
    detector_t det(resource);

    // context object
    typename detector_t::context ctx0{};

    brl_m_config brl_config{};
    edc_m_config edc_config{};

    if (n_edc_layers > edc_positions.size()) {
        throw std::invalid_argument(
            "ERROR: Too many endcap layers requested (max " +
            std::to_string(edc_positions.size()) + ")!");
    }
    if (n_brl_layers > brl_positions.size() - 1) {
        throw std::invalid_argument(
            "ERROR: Too many barrel layers requested (max " +
            std::to_string(brl_positions.size() - 1) + ")!");
    }
    // the radius of the endcaps and  the barrel section need to match
    if (n_edc_layers > 0 and
        std::fabs(brl_lay_sizes[n_brl_layers].second - edc_config.outer_r) >
            std::numeric_limits<scalar>::epsilon()) {
        throw std::invalid_argument(
            "ERROR: Barrel and endcap radii do not match!");
    }

    // beampipe
    dindex beampipe_idx = 0;
    add_beampipe(det, resource, ctx0, n_edc_layers, n_brl_layers, edc_lay_sizes,
                 brl_lay_sizes[0], brl_positions[0], brl_half_z,
                 edc_config.inner_r);

    if (n_edc_layers > 0) {
        edc_config.side = -1;
        // negative endcap layers
        add_endcap_detector<empty_vol_factory, edc_module_factory>(
            det, resource, ctx0, n_edc_layers, beampipe_idx, edc_lay_sizes,
            edc_positions, edc_config);

        // gap volume that connects barrel and neg. endcap
        dindex prev_vol_idx = det.volumes().back().index();
        prev_vol_idx = prev_vol_idx == 0 ? leaving_world : prev_vol_idx;
        dindex next_vol_idx = n_brl_layers == 0
                                  ? leaving_world
                                  : det.volumes().back().index() + 2;

        add_endcap_barrel_connection(
            det, resource, ctx0, edc_config.side, n_brl_layers, beampipe_idx,
            brl_lay_sizes, edc_config.inner_r, edc_config.outer_r,
            edc_lay_sizes[0].first, brl_half_z, next_vol_idx, prev_vol_idx);
    }
    if (n_brl_layers > 0) {
        // barrel
        add_barrel_detector<empty_vol_factory, brl_module_factory>(
            det, resource, ctx0, n_brl_layers, beampipe_idx, brl_half_z,
            brl_lay_sizes, brl_positions, brl_binning, brl_config);
    }
    if (n_edc_layers > 0) {
        // gap layer that connects barrel and pos. endcap
        edc_config.side = 1.;
        // innermost barrel layer volume id
        dindex prev_vol_idx =
            n_brl_layers == 0 ? leaving_world : 2 * n_edc_layers + 1;
        dindex next_vol_idx = prev_vol_idx == 1
                                  ? leaving_world
                                  : det.volumes().back().index() + 2;

        add_endcap_barrel_connection(
            det, resource, ctx0, edc_config.side, n_brl_layers, beampipe_idx,
            brl_lay_sizes, edc_config.inner_r, edc_config.outer_r, brl_half_z,
            edc_lay_sizes[0].first, prev_vol_idx, next_vol_idx);

        // positive endcap layers
        add_endcap_detector<empty_vol_factory, edc_module_factory>(
            det, resource, ctx0, n_edc_layers, beampipe_idx, edc_lay_sizes,
            edc_positions, edc_config);
    }

    return det;
}

}  // namespace detray