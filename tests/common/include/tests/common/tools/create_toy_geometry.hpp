/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/core/detector.hpp"
#include "tests/common/tools/detector_registry.hpp"

#include <vecmem/memory/host_memory_resource.hpp>

namespace detray {

/** Function that adds a cylinder portal.
 *
 * @param volume_id volume the portal should be added to
 * @param r raius of the cylinder
 * @param half_z z half length of the cylinder
 */
template <typename context_t, typename surface_container_t, typename mask_container_t, typename transform_container_t, typename edge_links, unsigned int cylinder_id = 3>
inline auto add_cylinder_portal(const dindex volume_id,
                         context_t& ctx,
                         surface_container_t& surfaces,
                         mask_container_t& masks,
                         transform_container_t& transforms,
                         const scalar r, const scalar lower_z, 
                         const scalar upper_z, const edge_links edge) {
    // translation
    point3 tsl{0., 0., 0};

    // add transform
    transforms[cylinder_id].emplace_back(ctx, tsl);

    // add mask
    masks.template add_mask<cylinder_id>(r, lower_z, upper_z);
    masks.template group<cylinder_id>().back().links() = edge;

    // add surface
    typename surface_container_t::value_type::value_type::mask_links mask_link{cylinder_id, masks.template size<cylinder_id>() - 1};
    surfaces[cylinder_id].emplace_back(transforms[cylinder_id].size(ctx) - 1, mask_link, volume_id, dindex_invalid);

    surfaces[cylinder_id].back().set_edge(edge);
}


/** Function that adds a disc portal.
 *
 * @param volume_id volume the portal should be added to
 * @param min_r lower radius of disc
 * @param max_r upper radius of disc
 * @param half_z z half length of the detector volume
 */
template <typename context_t, typename surface_container_t, typename mask_container_t, typename transform_container_t, typename edge_links, unsigned int disc_id = 4>
inline auto add_disc_portal(const dindex volume_id, 
                    context_t& ctx,
                    surface_container_t& surfaces, mask_container_t& masks,
                    transform_container_t& transforms, const scalar min_r,
                    const scalar max_r, const scalar z, 
                    const edge_links edge) {
    // translation
    point3 tsl{0., 0., z};

    // add transform
    transforms[disc_id].emplace_back(ctx, tsl);

    // add mask
    masks.template add_mask<disc_id>(min_r, max_r);
    masks.template group<disc_id>().back().links() = edge;

    // add surface
    typename surface_container_t::value_type::value_type::mask_links mask_link{disc_id, masks.template size<disc_id>() - 1};
    surfaces[disc_id].emplace_back(transforms[disc_id].size(ctx) - 1, mask_link,
        volume_id, dindex_invalid);

    surfaces[disc_id].back().set_edge(edge);
}

/** Creates a number of pixel modules for the a cylindrical barrel region.
 *
 * @tparam surface_t The surface type that contains the container indices
 * @tparam mask_t The geomterical boundaries of a surface (rectangles)
 *
 * @param m_half_x module half length in local x
 * @param m_half_y module half length in local y
 * @param m_tilt_phi phi tilt of the modules
 * @param layer_r radius at which the modules are positioned in the volume
 * @param radial_stagger module stagger in r
 * @param l_overlap the z overlap of modules next in z
 * @param binning phi, z bins of the surface grid
 *
 * @return a tuple that contains the surfaces (linking into the locally
 *         created container), the module transsforms and the surface masks.
 */
template <typename context_t, typename surface_container_t, typename mask_container_t, typename transform_container_t, unsigned int rectangle_id = 0>
inline auto create_barrel_modules(context_t& ctx, const dindex volume_id,
                           surface_container_t& surfaces,
                           mask_container_t& masks,
                           transform_container_t& transforms,
                           const scalar m_half_x = 8.4,
                           const scalar m_half_y = 36.,
                           const scalar m_tilt_phi = 0.14,
                           const scalar layer_r = 32.,
                           const scalar radial_stagger = 0.5,
                           const scalar l_overlap = 2.,
                           const std::pair<int, int> binning = {16, 14}) {
    /// mask index: type, range
    using mask_link_t = typename surface_container_t::value_type::value_type::mask_links;

    // Create the module centers

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
    scalar z_start =
        scalar{-0.5} * (n_z_bins - 1) * (scalar{2} * m_half_y - l_overlap);
    scalar z_step = scalar{2} * std::abs(z_start) / (n_z_bins - 1);

    // loop over the z bins
    for (size_t z_bin = 0; z_bin < size_t(n_z_bins); ++z_bin) {
        // prepare z and r
        scalar m_z = z_start + z_bin * z_step;
        scalar m_r = (z_bin % 2) != 0u ? layer_r - scalar{0.5} * radial_stagger
                                       : layer_r + scalar{0.5} * radial_stagger;
        for (size_t phiBin = 0; phiBin < size_t(n_phi_bins); ++phiBin) {
            // calculate the current phi value
            scalar m_phi = min_phi + phiBin * phi_step;
            m_centers.push_back(
                point3{m_r * std::cos(m_phi), m_r * std::sin(m_phi), m_z});
        }
    }

    // Create geometry data

    // First value is two for rectangle type, then index into local container
    mask_link_t m_id = {0, 0};

    for (auto& m_center : m_centers) {

        // Surfaces with the linking into the local containers
        m_id = {rectangle_id, masks.template size<rectangle_id>()};
        surfaces[rectangle_id].emplace_back(transforms[rectangle_id].size(ctx), m_id, volume_id, dindex_invalid);
        surfaces[rectangle_id].back().set_edge({volume_id, dindex_invalid});

        // The rectangle bounds for this module
        masks.template add_mask<rectangle_id>(m_half_x, m_half_y);
        masks.template group<rectangle_id>().back().links() = {volume_id, dindex_invalid};

        // Build the transform
        // The local phi
        scalar m_phi = algebra::getter::phi(m_center);
        // Local z axis is the normal vector
        vector3 m_local_z{std::cos(m_phi + m_tilt_phi),
                          std::sin(m_phi + m_tilt_phi), 0.};
        // Local x axis the normal to local y,z
        vector3 m_local_x{-std::sin(m_phi + m_tilt_phi),
                          std::cos(m_phi + m_tilt_phi), 0.};

        // Create the module transform
        transforms[rectangle_id].emplace_back(ctx, m_center, m_local_z, m_local_x);
    }
}

/** Helper method for positioning
  *
  * @param z is the z position of the ring
  * @param radius is the ring radius
  * @param phi_stagger is the radial staggering along phi
  * @param phi_sub_stagger is the overlap of the modules
  * @param n_phi_bins is the number of bins in phi
  *
  * @return a vector of the module positions in a ring
  */
inline auto module_positions_ring(scalar z,
                                  scalar radius,
                                  scalar phi_stagger,
                                  scalar phi_sub_stagger,
                                  int n_phi_bins) {
    // create and fill the positions
    std::vector<vector3> r_positions;
    r_positions.reserve(n_phi_bins);

    // prep work
    scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step = scalar{2} * pi / (n_phi_bins);
    double min_phi = -pi + 0.5 * phi_step;

    for (size_t iphi = 0; iphi < size_t(n_phi_bins); ++iphi) {
        // if we have a phi sub stagger presents
        double rzs = 0.;
        // phi stagger affects 0 vs 1, 2 vs 3 ... etc
        // -> only works if it is a %4
        // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
        if (phi_sub_stagger != 0. && !(n_phi_bins % 4)) {
            // switch sides
            if (!(iphi % 4)) {
                rzs = phi_sub_stagger;
            }
            else if (!((iphi + 1) % 4)) {
                rzs = -phi_sub_stagger;
            }
        }
        // the module phi
        double phi = min_phi + iphi * phi_step;
        // main z position depending on phi bin
        double rz = iphi % 2 ? z - 0.5 * phi_stagger : z + 0.5 * phi_stagger;
        r_positions.push_back(
            vector3{radius * std::cos(phi), radius * std::sin(phi), rz + rzs});
    }
    return r_positions;
}


/** Helper method for positioning rings in a disc
*
* @param z is the nominal z posiiton of the dis
* @param ring_stagger is the staggering of the different rings
* @param phi_stagger is the staggering on a ring in phi : it is even/odd
* @param phi_sub_stagger is the sub staggering on a ring in phi : it affects
* 0/4/8 and 3/6
* @param inner_r is the inner Radius for the disc
* @param outer_r is the outer Radius for the disc
* @param disc_binning is the binning setup in r (size of the vector), phi
* @param m_half_y is a pair of phibins and module length
* @param m_half_x_min_y The half lenght in X (at Y min) of the module
* @param m_half_x_max_y The half lenght in X (at Y max) of the module
* @param m_half_y The half lenght in Y of the module
* @param m_tilt The tilt out of the plane for discs
*
* @return vector of module positions of a ring
*/
template <typename context_t, typename surface_container_t, typename mask_container_t, typename transform_container_t, unsigned int trapezoid_id = 1>
void create_endcap_modules(context_t& ctx, const dindex volume_id,    
                           surface_container_t &surfaces,
                           mask_container_t &masks,
                           transform_container_t &transforms, scalar z,
                           const scalar ring_stagger = 0.0, 
                           const std::vector<scalar> phi_stagger = {4.0, 4.0},
                           const std::vector<scalar> phi_sub_stagger = {0., 0.},
                           const scalar inner_r = 27., const scalar outer_r = 180.,
                           const std::vector<size_t>& disc_binning = {40, 68},
                           const std::vector<scalar>& m_half_y = {36., 36.},
                           const std::vector<scalar>& m_half_x_min_y = {8.4, 8.4}, 
                           const std::vector<scalar>& m_half_x_max_y = {12.4, 12.4}, 
                           const std::vector<scalar>& m_tilt = {0., 0.},
                           int side = 1) {
    using mask_link_t = typename surface_container_t::value_type::value_type::mask_links;
    // calculate the radii of the rings
    std::vector<scalar> radii;
    // calculate the radial borders
    //std::vector<scalar> radial_boarders;
    // the radial span of the disc
    scalar delta_r = outer_r - inner_r;

    // Only one ring
    if (disc_binning.size() == 1) {
        radii.push_back(scalar{0.5} * (inner_r + outer_r));
        //radial_boarders = {inner_r, outer_r};
    }
    else {
        // sum up the total length of the modules along r
        scalar tot_length = 0;
        for (auto& m_hlength : m_half_y) {
            tot_length += scalar{2} * m_hlength + 0.5;
        }
        // now calculate the overlap (equal pay)
        scalar r_overlap = (tot_length - delta_r) / (m_half_y.size() - 1);
        // and now fill the radii and gaps
        scalar prev_r = inner_r;
        scalar prev_hl = 0.;
        scalar prev_ol = 0.;
        // remember the radial boarders
        //radial_boarders.push_back(inner_r);
        for (auto& m_hlength : m_half_y) {
            // calculate the radius
            radii.push_back(prev_r + prev_hl - prev_ol + m_hlength);
            prev_r = radii.back();
            prev_ol = r_overlap;
            prev_hl = m_hlength;
            // and register the radial boarder
            //radial_boarders.push_back(prev_r + scalar{2} * prev_hl - scalar{0.5} * prev_ol);
        }
    }

    // now build the modules in every ring
    for (size_t ir = 0; ir < radii.size(); ++ir) {
        // generate the z value
        // convention inner ring is closer to origin : makes sense
        double rz = radii.size() == 1
                        ? z
                        : (ir % 2 ? z + scalar{0.5} * ring_stagger : z - scalar{0.5} * ring_stagger);
        // fill the ring module positions
        double ps_stagger = phi_sub_stagger.size() ? phi_sub_stagger[ir] : 0.;

        std::vector<point3> r_postitions = module_positions_ring(rz, radii[ir], 
                                                    phi_stagger[ir],
                                                    ps_stagger,
                                                    disc_binning[ir]);

        // Build the geometrical objects
        for (const auto &m_position : r_postitions) {
            // trapezoid mask
            mask_link_t mask_link{trapezoid_id, masks.template size<trapezoid_id>()};
            masks.template add_mask<trapezoid_id>(m_half_x_min_y[ir], m_half_x_max_y[ir], m_half_y[ir]);
            // The links will be updated to the volume later
            masks.template group<trapezoid_id>().back().links() = {volume_id, dindex_invalid};

            // Surfaces with the linking into the local containers
            surfaces[trapezoid_id].emplace_back(transforms[trapezoid_id].size(ctx), mask_link, volume_id, dindex_invalid);
            surfaces[trapezoid_id].back().set_edge({volume_id, dindex_invalid});

            // the module transform from the position
            double m_phi = algebra::getter::phi(m_position);
            // the center position of the modules
            point3 m_center{static_cast<scalar>(side) * m_position};
            // the rotation matrix of the module
            vector3 m_local_y{std::cos(m_phi), std::sin(m_phi), 0.};
            // take different axis to have the same readout direction
            vector3 m_local_z{0., 0., side * 1.};
            vector3 m_local_x = algebra::vector::cross(m_local_y, m_local_z);

            // Create the module transform
            transforms[trapezoid_id].emplace_back(ctx, m_center, m_local_z, m_local_x);
        }
    }
}


/** Helper method for positioning rings in a disc
*
* @param z is the nominal z posiiton of the dis
* @param ring_stagger is the staggering of the different rings
* @param phi_stagger is the staggering on a ring in phi : it is even/odd
* @param phi_sub_stagger is the sub staggering on a ring in phi : it affects
* 0/4/8 and 3/6
* @param inner_r is the inner Radius for the disc
* @param outer_r is the outer Radius for the disc
* @param disc_binning is the binning setup in r (size of the vector), phi
* @param m_half_y is a pair of phibins and module length
* @param m_half_x_min_y The half lenght in X (at Y min) of the module
* @param m_half_x_max_y The half lenght in X (at Y max) of the module
* @param m_half_y The half lenght in Y of the module
* @param m_tilt The tilt out of the plane for discs
*
* @return vector of module positions of a ring
*/
template <typename detector_t, typename context_t>
inline auto add_beampipe(detector_t &det, vecmem::memory_resource &resource,context_t &ctx, const std::vector<std::pair<scalar, scalar>> &edc_vol_sizes, const std::pair<scalar, scalar> &beampipe_vol_size, const scalar beampipe_r, const scalar brl_half_z, const scalar edc_inner_r) {

    auto constexpr cylinder_id = detector_t::e_cylinder3;

    const dindex inv_sf_finder = dindex_invalid;
    const dindex leaving_world = dindex_invalid;

    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_container transforms = {resource};

    auto& beampipe = det.new_volume({beampipe_vol_size.first, beampipe_vol_size.second, -edc_vol_sizes[2].second, edc_vol_sizes[2].second, -M_PI, M_PI});
    const auto beampipe_idx = beampipe.index();

    // This is the beampipe surface
    // identity
    transforms[cylinder_id].emplace_back(ctx);
    // add mask
    masks.template add_mask<cylinder_id>(beampipe_r, 
                                         -edc_vol_sizes[2].second, 
                                         edc_vol_sizes[2].second);
    masks.template group<cylinder_id>().back().links() = {beampipe_idx, inv_sf_finder};
    // add surface
    typename detector_t::mask_link mask_link = {cylinder_id, masks.template size<cylinder_id>() - 1};
    surfaces[cylinder_id].emplace_back(
        transforms[cylinder_id].size(ctx) - 1, std::move(mask_link),
        beampipe_idx, dindex_invalid);
    surfaces[cylinder_id].back().set_edge({beampipe_idx, inv_sf_finder});

    // negative and positive, outer portal surfaces
    // cylinder portals for all volumes
    //negative endcap
    typename detector_t::edge_type edge = {beampipe_idx + 1, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, -edc_vol_sizes[2].second, -edc_vol_sizes[2].first, edge);
    edge = {beampipe_idx + 2, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, -edc_vol_sizes[2].first, -edc_vol_sizes[1].second, edge);
    edge = {beampipe_idx + 3, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, -edc_vol_sizes[1].second, -edc_vol_sizes[1].first, edge);
    edge = {beampipe_idx + 4, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, -edc_vol_sizes[1].first, -edc_vol_sizes[0].second, edge);
    edge = {beampipe_idx + 5, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, -edc_vol_sizes[0].second, -edc_vol_sizes[0].first, edge);
    edge = {beampipe_idx + 6, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, -edc_vol_sizes[0].first, -brl_half_z, edge);
    // barrel
    edge = {beampipe_idx + 7, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, -brl_half_z, brl_half_z, edge);
    // positive endcap
    edge = {beampipe_idx + 14, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, brl_half_z, edc_vol_sizes[0].first, edge);
    edge = {beampipe_idx + 15, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, edc_vol_sizes[0].first, edc_vol_sizes[0].second, edge);
    edge = {beampipe_idx + 16, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, edc_vol_sizes[0].second, edc_vol_sizes[1].first, edge);
    edge = {beampipe_idx + 17, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, edc_vol_sizes[1].first, edc_vol_sizes[1].second, edge);
    edge = {beampipe_idx + 18, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, edc_vol_sizes[1].second, edc_vol_sizes[2].first, edge);
    edge = {beampipe_idx + 19, inv_sf_finder};
    add_cylinder_portal(beampipe_idx, ctx, surfaces, masks, transforms, edc_inner_r, edc_vol_sizes[2].first, edc_vol_sizes[2].second, edge);
    // discs
    edge = {leaving_world, inv_sf_finder};
    add_disc_portal(beampipe_idx, ctx, surfaces, masks, transforms, beampipe_vol_size.first, beampipe_vol_size.second, -edc_vol_sizes[2].second, edge);
    edge = {leaving_world, inv_sf_finder};
    add_disc_portal(beampipe_idx, ctx, surfaces, masks, transforms, beampipe_vol_size.first, beampipe_vol_size.second, edc_vol_sizes[2].second, edge);

    det.add_objects(ctx, beampipe, surfaces, masks, transforms);
}


/// @return an endcap volume
template<typename detector_t>
void add_endcap_volume(detector_t &det, vecmem::memory_resource &resource,
                       typename detector_t::context &ctx, const scalar side, const scalar lay_neg_z, const scalar lay_pos_z,
                       const scalar edc_inner_r, const scalar edc_outer_r, 
                       const std::vector<typename detector_t::edge_type> &edges, const bool is_gap = false, 
                       const scalar edc_position = 600.,
                       const scalar ring_stagger = 0.0, 
                       const std::vector<scalar> m_phi_stagger = {4.0, 4.0},
                       const std::vector<scalar> m_phi_sub_stagger = {0.5, 0.},
                       const std::vector<size_t>& disc_binning = {40, 68},
                       const std::vector<scalar>& m_half_y = {36., 36.},
                       const std::vector<scalar> m_half_x_min_y = {8.4, 8.4}, 
                       const std::vector<scalar> m_half_x_max_y = {12.4, 12.4}, 
                       const std::vector<scalar> m_tilt = {0., 0.}) {
    // build the outermost volume
    const scalar edc_lower_z = std::min(side*lay_neg_z, side*lay_pos_z);
    const scalar edc_upper_z = std::max(side*lay_neg_z, side*lay_pos_z);;

    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_container transforms = {resource};

    auto& edc_volume = det.new_volume({edc_inner_r, edc_outer_r, edc_lower_z, edc_upper_z, -M_PI, M_PI});

    if (not is_gap) {
        // create disc module surfaces
        create_endcap_modules(ctx, edc_volume.index(), 
                              surfaces, masks, transforms,
                              side * edc_position, ring_stagger, m_phi_stagger,
                              m_phi_sub_stagger, edc_inner_r, edc_outer_r,
                              disc_binning, m_half_y, m_half_x_min_y,
                              m_half_x_max_y, m_tilt);
    }

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(edc_volume.index(), ctx, surfaces, masks, transforms, edc_inner_r, edc_lower_z, edc_upper_z, edges[0]);
    add_cylinder_portal(edc_volume.index(), ctx, surfaces, masks, transforms, edc_outer_r, edc_lower_z, edc_upper_z, edges[1]);
    add_disc_portal(edc_volume.index(), ctx, surfaces, masks, transforms, edc_inner_r, edc_outer_r, edc_lower_z, edges[2]);
    add_disc_portal(edc_volume.index(), ctx, surfaces, masks, transforms, edc_inner_r, edc_outer_r, edc_upper_z, edges[3]);

    det.add_objects(ctx, edc_volume, surfaces, masks, transforms);
}

template<typename detector_t>
void add_barrel_volume(detector_t &det, vecmem::memory_resource &resource,
                       typename detector_t::context &ctx, const scalar lay_inner_r, const scalar lay_outer_r,
                       const scalar brl_half_z, 
                       const std::vector<typename detector_t::edge_type> &edges, const bool is_gap = false, 
                       const scalar m_half_x = 8.4,
                       const scalar m_half_y = 36.,
                       const scalar m_tilt_phi = 0.14,
                       const scalar layer_r = 32.,
                       const scalar m_radial_stagger = 0.5,
                       const scalar m_long_overlap = 2.,
                       const std::pair<int, int> m_binning = {16, 14}) {

    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_container transforms = {resource};

    auto& brl_volume = det.new_volume({lay_inner_r, lay_outer_r, -brl_half_z, brl_half_z, -M_PI, M_PI});

    if (not is_gap) {
        // create disc module surfacesinline auto 
        create_barrel_modules(ctx, brl_volume.index(),
                              surfaces, masks, transforms, m_half_x, 
                              m_half_y, m_tilt_phi, layer_r,
                              m_radial_stagger, m_long_overlap, m_binning);
    }

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_volume.index(), ctx, surfaces, masks, transforms, lay_inner_r, -brl_half_z, brl_half_z, edges[0]);
    add_cylinder_portal(brl_volume.index(), ctx, surfaces, masks, transforms, lay_outer_r, -brl_half_z, brl_half_z, edges[1]);
    add_disc_portal(brl_volume.index(), ctx, surfaces, masks, transforms, lay_inner_r, lay_outer_r, -brl_half_z, edges[2]);
    add_disc_portal(brl_volume.index(), ctx, surfaces, masks, transforms, lay_inner_r, lay_outer_r, brl_half_z, edges[3]);

    det.add_objects(ctx, brl_volume, surfaces, masks, transforms);
}


/// @return an endcap volume
/*std::shared_ptr<DetectorVolume> createEndcapVolume(
    scalar volumeMinR, scalar volumeMaxR, scalar volumeMinZ,
    scalar volumeMaxZ, const std::string& volumeName = "SingleLayerVolume",
    scalar moduleHalfXminY = 8.4, scalar moudleHalfXmaxY = 12.4,
    scalar moduleHalfY = 32., scalar moduleTilt = 0.,
    scalar ringRadius = 40., scalar zStagger = 2, int nPhi = 40) {
  // Place the ring into the middle
  scalar ringZ = 0.5 * (volumeMinZ + volumeMaxZ);

  auto volumeSurfaces =
      surfacesRing(moduleHalfXminY, moudleHalfXmaxY, moduleHalfY, moduleTilt,
                   ringRadius, ringZ, zStagger, nPhi);

  // Create the volume bounds
  auto volumeBounds = std::make_unique<CylinderVolumeBounds>(
      volumeMinR, volumeMaxR, 0.5 * (volumeMaxZ - volumeMinZ));

  SurfaceLinks volumeSurfaceLinks = AllSurfaces{};

  std::vector<SurfaceLinks> portalSurfaceLinks = {AllSurfaces{}, AllSurfaces{},
                                                  AllSurfaces{}};

  if (volumeMinR > 0.) {
    portalSurfaceLinks.push_back(AllSurfaces{});
  }

  auto volumeTransform = Transform3::Identity();
  volumeTransform.pretranslate(Vector3(0., 0., ringZ));

  return DetectorVolume::makeShared(
      volumeTransform, std::move(volumeBounds), std::move(volumeSurfaces),
      std::move(volumeSurfaceLinks), std::move(portalSurfaceLinks), volumeName);
}*/


/// Helper method to create a central detector
///
/// Keep track of the central half length with
/// @param detectorRmin inner radius of detector
/// @param detectorRmax outer radius of detector
/// @param zToCentral is the distance to central to
/// @param side is the side of the endcap detector
/// @param detectorName is the detector name prescript
///
/// @return a central detector volume
/*std::shared_ptr<DetectorVolume> createEndcapDetector(
    scalar detectorRmin = 0., scalar detectorRmax = 80.,
    scalar zToCentral = 500., int side = 1) {
  scalar layerThickness = 5.;
  scalar gapThickness = 50.;

  std::string firstLayerName = (side > 0) ? "Layer0" : "Layer1";
  std::string gapName = "Gap";
  std::string secondLayerName = (side > 0) ? "Layer1" : "Layer0";
  std::string sideTag = (side > 0) ? "Pos" : "Neg";

  // Place the first layer
  scalar oneZ = side * (zToCentral);
  scalar twoZ = side * (zToCentral + layerThickness);
  auto firstLayer = createEndcapVolume(
      detectorRmin, detectorRmax, std::min(oneZ, twoZ), std::max(oneZ, twoZ),
      detectorName + firstLayerName + sideTag);

  // Adapt for the gap & build
  oneZ = side * (zToCentral + layerThickness);
  twoZ = side * (zToCentral + layerThickness + gapThickness);
  Transform3 gapTransform = Transform3::Identity();
  gapTransform.pretranslate(Vector3(0., 0., 0.5 * (oneZ + twoZ)));
  auto gapBounds = std::make_unique<CylinderVolumeBounds>(
      detectorRmin, detectorRmax, std::abs(0.5 * (oneZ - twoZ)));
  auto gap = DetectorVolume::makeShared(gapTransform, std::move(gapBounds),
                                        detectorName + gapName + sideTag);

  // Adapt for the second layer
  oneZ = side * (zToCentral + layerThickness + gapThickness);
  twoZ = side * (zToCentral + 2 * layerThickness + gapThickness);
  auto secondLayer = createEndcapVolume(
      detectorRmin, detectorRmax, std::min(oneZ, twoZ), std::max(oneZ, twoZ),
      detectorName + secondLayerName + sideTag);

  std::vector<std::shared_ptr<DetectorVolume>> endcapVolumes;
  if (side > 0) {
    endcapVolumes = {firstLayer, gap, secondLayer};
  } else {
    endcapVolumes = {secondLayer, gap, firstLayer};
  }
  // Container in Z
  return CylindricalContainerHelper::containerInZ(
      std::move(endcapVolumes),
      detectorName + std::string("TwoLayers") + sideTag);
}*/

/** Builds a simple detray geometry of the innermost tml layers. It contains:
 *
 * - a beampipe (r = 19mm, half_z = 825mm)
 * - a first layer (r_min = 27mm, r_max = 38mm, half_z = 500mm) with 224
 *   rectangular (half_x = 8.4mm, half_y = 36mm) modules at r = 32mm
 * - an empty layer (r_min = 38mm, r_max = 64mm, half_z = 500mm)
 * - a second layer (r_min = 64mm, r_max = 80mm, half_z = 500mm) with 448
 *   rectangular (half_x = 8.4mm, half_y = 36mm) modules at r = 72mm.
 *
 * @returns a tuple containing the geometry objects collections: [volumes,
 *          surfaces, transforms, disc masks (neg/pos portals), cylinder masks
 *          (inner/outer portals), rectangle masks (modules)]
 */
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
    using volume = typename detector_t::volume_type;
    using surface = typename detector_t::surface_type;
    using mask_container = typename detector_t::mask_container;
    using transform_store = typename detector_t::transform_store;
    using transform_container = typename detector_t::transform_container;
    using surface_container = typename detector_t::surface_filling_container;

    /** source link */
    const dindex inv_sf_finder = dindex_invalid;
    /** Leaving world */
    const dindex leaving_world = dindex_invalid;

    //
    // barrel
    //
    const scalar brl_half_z = 500.;
    const std::vector<scalar> brl_positions = {19., 32., 72., 116., 172.};
    const std::vector<std::pair<scalar, scalar>> brl_vol_sizes = {{0., 27.}, {27., 38.}, {64., 80.}, {108., 124.}, {164., 180.}};
    const scalar brl_radial_stagger = 0.5;//2.;
    const scalar brl_l_overlap = 2.;//5.;
    const std::vector<std::pair<int, int>> brl_binning = {{0., 0.}, {16, 14}, {32, 14}, {52, 14}, {78, 14}};
    // module parameters
    const scalar brl_half_x = 8.4;
    const scalar brl_half_y = 36.;
    const scalar brl_tilt_phi = 0.14;//0.145;

    //
    // endcaps
    //
    const std::vector<scalar> edc_positions = {600., 700., 820., 960., 1100., 1300., 1500.};
    const std::vector<std::pair<scalar, scalar>> edc_vol_sizes = {{595., 605.}, {695., 705.}, {815., 825.}, {955., 965.}, {1095., 1105.}, {1295., 1305.}, {1495., 1505.}};
    const scalar edc_ring_stagger = 1.0;
    // Parameters for both rings of modules
    const std::vector<scalar> edc_phi_stagger = {4.0, 4.0};
    const std::vector<scalar> edc_phi_sub_stagger = {0.5, 0.};
    const scalar edc_inner_r = 27.;
    const scalar edc_outer_r = 180.;
    const std::vector<size_t>& edc_disc_binning = {40, 68};
    // module params
    const std::vector<scalar>& edc_half_y = {36., 36.};
    const std::vector<scalar> edc_half_x_min_y = {8.4, 8.4};
    const std::vector<scalar> edc_half_x_max_y = {12.4, 12.4};
    const std::vector<scalar> edc_tilt = {0., 0.};

    // create detector
    detector_t det(resource);

    // context objects
    typename transform_store::context ctx0;

    add_beampipe(det, resource, ctx0, edc_vol_sizes, brl_vol_sizes[0], brl_positions[0], brl_half_z, edc_inner_r);

    dindex beampipe_idx = 0;
    //
    // negative endcap
    //

    //
    // first layer: Contains modules
    //
    int side = -1;
    bool is_gap = true;

    dindex next_vol_idx = beampipe_idx + 2;
    std::vector<edge_links> edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {leaving_world, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[2].second, edc_vol_sizes[2].first, edc_inner_r, edc_outer_r, edges, !is_gap, edc_positions[2], edc_ring_stagger, edc_phi_stagger, edc_phi_sub_stagger, edc_disc_binning, edc_half_y, edc_half_x_min_y, edc_half_x_max_y, edc_tilt);

    //
    // gap: Connect portals (first layer, second layer)
    //
    dindex prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[2].first, edc_vol_sizes[1].second, edc_inner_r, edc_outer_r, edges, is_gap);

    //
    // second layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[1].second, edc_vol_sizes[1].first, edc_inner_r, edc_outer_r, edges, !is_gap, edc_positions[1], edc_ring_stagger, edc_phi_stagger, edc_phi_sub_stagger, edc_disc_binning, edc_half_y, edc_half_x_min_y, edc_half_x_max_y, edc_tilt);

    //
    // gap layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[1].first, edc_vol_sizes[0].second, edc_inner_r, edc_outer_r, edges, is_gap);

    //
    // third layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[0].second, edc_vol_sizes[0].first, edc_inner_r, edc_outer_r, edges, !is_gap, edc_positions[0], edc_ring_stagger, edc_phi_stagger, edc_phi_sub_stagger, edc_disc_binning, edc_half_y, edc_half_x_min_y, edc_half_x_max_y, edc_tilt);

    //
    // final gap layer connects to barrel
    //
    typename detector_t::surface_filling_container surfaces = {};
    typename detector_t::mask_container masks = {resource};
    typename detector_t::transform_container transforms = {resource};

    scalar gap_neg_z = side * edc_vol_sizes[0].first;
    scalar gap_pos_z = side * brl_half_z;

    auto& final_gap = det.new_volume({edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z, -M_PI, M_PI});

    dindex final_gap_idx = det.volumes().back().index();
    prev_vol_idx = final_gap_idx - 1;
    next_vol_idx = prev_vol_idx + 2;

    typename detector_t::edge_type edge = {beampipe_idx, inv_sf_finder};
    add_cylinder_portal(final_gap_idx, ctx0, surfaces, masks, transforms,edc_inner_r, gap_neg_z, gap_pos_z, edge);
    edge = {leaving_world, inv_sf_finder};
    add_cylinder_portal(final_gap_idx, ctx0, surfaces, masks, transforms,edc_outer_r, gap_neg_z, gap_pos_z, edge);
    edge = {prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,edc_inner_r, edc_outer_r, gap_neg_z, edge);
    edge = {next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,brl_vol_sizes[1].first, brl_vol_sizes[1].second, gap_pos_z, edge);
    edge = {++next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,brl_vol_sizes[1].second, brl_vol_sizes[2].first, gap_pos_z, edge);
    edge = {++next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,brl_vol_sizes[2].first, brl_vol_sizes[2].second, gap_pos_z, edge);
    edge = {++next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,brl_vol_sizes[2].second, brl_vol_sizes[3].first, gap_pos_z, edge);
    edge = {++next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,brl_vol_sizes[3].first, brl_vol_sizes[3].second, gap_pos_z, edge);
    edge = {++next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,brl_vol_sizes[3].second, brl_vol_sizes[4].first, gap_pos_z, edge);
    edge = {++next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, surfaces, masks, transforms,brl_vol_sizes[4].first, brl_vol_sizes[4].second, gap_pos_z, edge);

    det.add_objects(ctx0, final_gap, surfaces, masks, transforms);

    //
    // barrel
    //

    //
    // first layer
    //
    scalar brl_inner_r = brl_vol_sizes[1].first;
    scalar brl_outer_r = brl_vol_sizes[1].second;

    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}, {6, inv_sf_finder}, {14, inv_sf_finder}};

    add_barrel_volume(det, resource, ctx0, brl_inner_r, brl_outer_r, brl_half_z, edges, !is_gap, brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[1], brl_radial_stagger, brl_l_overlap, brl_binning[1]);

    //
    // gap layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}, {6, inv_sf_finder}, {14, inv_sf_finder}};

    add_barrel_volume(det, resource, ctx0, brl_vol_sizes[1].second, brl_vol_sizes[2].first, brl_half_z, edges, is_gap);

    //
    // second layer
    //
    brl_inner_r = brl_vol_sizes[2].first;
    brl_outer_r = brl_vol_sizes[2].second;

    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}, {6, inv_sf_finder}, {14, inv_sf_finder}};

    add_barrel_volume(det, resource, ctx0, brl_inner_r, brl_outer_r, brl_half_z, edges, !is_gap, brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[2], brl_radial_stagger, brl_l_overlap, brl_binning[2]);

    //
    // gap layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}, {6, inv_sf_finder}, {14, inv_sf_finder}};

    add_barrel_volume(det, resource, ctx0, brl_vol_sizes[2].second, brl_vol_sizes[3].first, brl_half_z, edges, is_gap);

    //
    // third layer
    //
    brl_inner_r = brl_vol_sizes[3].first;
    brl_outer_r = brl_vol_sizes[3].second;

    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}, {6, inv_sf_finder}, {14, inv_sf_finder}};

    add_barrel_volume(det, resource, ctx0, brl_inner_r, brl_outer_r, brl_half_z, edges, !is_gap, brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[3], brl_radial_stagger, brl_l_overlap, brl_binning[3]);

    //
    // gap layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}, {6, inv_sf_finder}, {14, inv_sf_finder}};

    add_barrel_volume(det, resource, ctx0, brl_vol_sizes[3].second, brl_vol_sizes[4].first, brl_half_z, edges, is_gap);

    //
    // fourth layer
    //
    brl_inner_r = brl_vol_sizes[4].first;
    brl_outer_r = brl_vol_sizes[4].second;

    prev_vol_idx = det.volumes().back().index();
    edges = {{prev_vol_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {6, inv_sf_finder}, {14, inv_sf_finder}};

    add_barrel_volume(det, resource, ctx0, brl_inner_r, brl_outer_r, brl_half_z, edges, !is_gap, brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[4], brl_radial_stagger, brl_l_overlap, brl_binning[4]);

    //
    // positive endcap
    //

    side = 1.;

    //
    // final gap layer connects to barrel
    //
    typename detector_t::surface_filling_container pos_edc_surfaces = {};
    typename detector_t::mask_container pos_edc_masks = {resource};
    typename detector_t::transform_container pos_edc_transforms = {resource};

    gap_neg_z = side * brl_half_z;
    gap_pos_z = side * edc_vol_sizes[0].first;

    auto& pos_final_gap = det.new_volume({edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z, -M_PI, M_PI});

    final_gap_idx = det.volumes().back().index();
    prev_vol_idx = 7;
    next_vol_idx = final_gap_idx + 1;

    edge = {beampipe_idx, inv_sf_finder};
    add_cylinder_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, edc_inner_r, gap_neg_z, gap_pos_z, edge);
    edge = {leaving_world, inv_sf_finder};
    add_cylinder_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, edc_outer_r, gap_neg_z, gap_pos_z, edge);
    edge = {prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, brl_vol_sizes[1].first, brl_vol_sizes[1].second, gap_neg_z, edge);
    edge = {++prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, brl_vol_sizes[1].second, brl_vol_sizes[2].first, gap_neg_z, edge);
    edge = {++prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, brl_vol_sizes[2].first, brl_vol_sizes[2].second, gap_neg_z, edge);
    edge = {++prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, brl_vol_sizes[2].second, brl_vol_sizes[3].first, gap_neg_z, edge);
    edge = {++prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, brl_vol_sizes[3].first, brl_vol_sizes[3].second, gap_neg_z, edge);
    edge = {++prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, brl_vol_sizes[3].second, brl_vol_sizes[4].first, gap_neg_z, edge);
    edge = {++prev_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, brl_vol_sizes[4].first, brl_vol_sizes[4].second, gap_neg_z, edge);
    edge = {next_vol_idx, inv_sf_finder};
    add_disc_portal(final_gap_idx, ctx0, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms, edc_inner_r, edc_outer_r, gap_pos_z, edge);

    det.add_objects(ctx0, pos_final_gap, pos_edc_surfaces, pos_edc_masks, pos_edc_transforms);

    //
    // first layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[0].second, edc_vol_sizes[0].first, edc_inner_r, edc_outer_r, edges, !is_gap, edc_positions[0], edc_ring_stagger, edc_phi_stagger, edc_phi_sub_stagger, edc_disc_binning, edc_half_y, edc_half_x_min_y, edc_half_x_max_y, edc_tilt);

    //
    // gap layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[1].first, edc_vol_sizes[0].second, edc_inner_r, edc_outer_r, edges, is_gap);

    //
    // second layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[1].second, edc_vol_sizes[1].first, edc_inner_r, edc_outer_r, edges, !is_gap, edc_positions[1], edc_ring_stagger, edc_phi_stagger, edc_phi_sub_stagger, edc_disc_binning, edc_half_y, edc_half_x_min_y, edc_half_x_max_y, edc_tilt);

    //
    // gap layer
    //
    prev_vol_idx = det.volumes().back().index();
    next_vol_idx = prev_vol_idx + 2;
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {next_vol_idx, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[2].first, edc_vol_sizes[1].second, edc_inner_r, edc_outer_r, edges, is_gap);

    //
    // third layer
    //
    prev_vol_idx = det.volumes().back().index();
    edges = {{beampipe_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}, {prev_vol_idx, inv_sf_finder}, {leaving_world, inv_sf_finder}};

    add_endcap_volume(det, resource, ctx0, side, edc_vol_sizes[2].first, edc_vol_sizes[2].second, edc_inner_r, edc_outer_r, edges, !is_gap, edc_positions[2], edc_ring_stagger, edc_phi_stagger, edc_phi_sub_stagger, edc_disc_binning, edc_half_y, edc_half_x_min_y, edc_half_x_max_y, edc_tilt);

    return std::move(det);
}

}  // namespace detray