/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <cmath>
#include <map>
#include <string>

#include "detray/core/detector.hpp"
#include "detray/geometry/surface_base.hpp"
#include "detray/geometry/volume.hpp"
#include "detray/io/csv_io.hpp"
#include "detray/masks/masks.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

/** Keeps the relevant csv file names */
struct detector_input_files {
    std::string det_name, surface, layer_volume, surface_grid,
        surface_grid_entries;
};

/// open data detector file names
detector_input_files odd_files = {"odd", "odd.csv", "odd-layer-volumes.csv",
                                  "odd-surface-grids.csv", ""};

/// track ml detector file names
detector_input_files tml_files = {"tml", "tml.csv", "tml-layer-volumes.csv",
                                  "tml-surface-grids.csv", ""};

/** Read a detector from csv files */
auto read_from_csv(detector_input_files& files,
                   vecmem::memory_resource& resource) {
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure(
            "Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string surfaces = data_directory + files.surface;
    std::string volumes = data_directory + files.layer_volume;
    std::string grids = data_directory + files.surface_grid;
    std::string grid_entries = files.surface_grid_entries;
    std::map<detray::dindex, std::string> name_map{};

    auto d =
        detray::detector_from_csv<>(files.det_name, surfaces, volumes, grids,
                                    grid_entries, name_map, resource);

    return std::make_pair<decltype(d), decltype(name_map)>(std::move(d),
                                                           std::move(name_map));
};

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
template <typename surface_t, typename mask_t>
inline auto create_modules(scalar m_half_x = 8.4, scalar m_half_y = 36.,
                           scalar m_tilt_phi = 0.145, scalar layer_r = 32.,
                           scalar radial_stagger = 2., scalar l_overlap = 5.,
                           const std::pair<int, int> binning = {16, 14}) {

    /// mask index: type, range
    using mask_index = detray::darray<detray::dindex, 2>;

    // Prepare the return container
    using surface_container = detray::dvector<surface_t>;
    using trfs_container = detray::dvector<transform3>;
    using mask_container = detray::dvector<mask_t>;

    surface_container surfaces;
    trfs_container transforms;
    mask_container masks;

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

    // loop over the bins
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
    mask_index m_id = {2, 0};

    for (auto& m_center : m_centers) {

        // Surfaces with the linking into the local containers
        m_id = {2, masks.size()};
        surfaces.emplace_back(transforms.size(), m_id, detray::dindex_invalid,
                              detray::dindex_invalid);

        // The rectangle bounds for this module
        masks.emplace_back(m_half_x, m_half_y);

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
        transforms.emplace_back(m_center, m_local_z, m_local_x);
    }

    return std::make_tuple<surface_container, trfs_container, mask_container>(
        std::move(surfaces), std::move(transforms), std::move(masks));
}

// Types for toy geometry
struct object_registry {
    // Known primitives
    enum id : unsigned int {
        e_object_types = 1,
        e_surface = 0,
        e_any = 0,
        e_unknown = 2,
    };

    template <typename value_type = void>
    static constexpr auto get() {
        return e_surface;
    }
};

/** Builds a simple detray geometry of the innermost tml layers. It contains:
 *
 * - a beampipe (r = 27mm, half_z = 500mm)
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
auto toy_geometry() {

    // Volume type
    using volume_type =
        detray::volume<object_registry, dindex_range, detray::darray>;
    /// volume index: volume the surface belongs to
    using volume_index = detray::dindex;
    /// transform link: transform entry belonging to surface
    using transform_link = detray::dindex;
    /// mask index: type, range
    using mask_index = detray::darray<detray::dindex, 2>;
    /// volume links: next volume, next (local) object finder
    using edge_links = detray::darray<detray::dindex, 2>;
    /// source link
    using source_link = detray::dindex;

    // We have disc, cylinder and rectangle surfaces
    using disc = detray::ring2<detray::planar_intersector, __plugin::cartesian2,
                               short, 0>;

    using cylinder = detray::cylinder3<false, detray::cylinder_intersector,
                                       __plugin::cylindrical2, short, 1>;

    using rectangle = detray::rectangle2<detray::planar_intersector,
                                         __plugin::cartesian2, short, 2>;

    // The surface type, both for volume portals and contained detector
    // surfaces
    using surface = detray::surface_base<transform_link, mask_index,
                                         volume_index, source_link, edge_links>;

    // The geometry data
    using volume_container = detray::dvector<volume_type>;
    using surface_container = detray::dvector<surface>;
    using transf_container = detray::dvector<transform3>;
    using disc_container = detray::dvector<disc>;
    using cylinder_container = detray::dvector<cylinder>;
    using rectangle_container = detray::dvector<rectangle>;

    /** Volumes */
    volume_container volumes = {};
    /** Surfaces, including portals */
    surface_container surfaces = {};
    /** Surface transforms */
    transf_container transforms = {};
    /** Disc masks for neg/pos volume boundaries (portals) */
    disc_container discs = {};
    /** Cylinder masks for inner/outer volume boundaries (portals) */
    cylinder_container cylinders = {};
    /** Rectangle masks for detector surfaces */
    rectangle_container rectangles = {};
    // mask index for surfaces
    mask_index m_id = {};
    /** source link */
    const dindex inv_sf_finder = dindex_invalid;
    ;
    /** Leaving world */
    const dindex leaving_world = dindex_invalid;

    // parameters
    scalar detector_half_z = 500.;
    scalar beampipe_r = 27.;
    scalar first_layer_outer_r = 38.;
    scalar second_layer_inner_r = 64.;
    scalar second_layer_outer_r = 80.;

    /** Add a single barrel layer volume to an existing collection.
     *
     * @param min_r minimal radius of volume
     * @param max_r maximal radius of volume
     * @param half_z half length in z of volume
     *
     * @return a detray cylinder volume
     */
    auto add_cylinder_volume = [&](scalar min_r, scalar max_r,
                                   scalar half_z) -> volume_type& {
        // The volume bounds
        detray::darray<scalar, 6> bounds = {min_r,  max_r, -half_z,
                                            half_z, -M_PI, M_PI};
        // Add the new volume to the collection
        auto& vol = volumes.emplace_back(bounds);
        vol.set_index(volumes.size() - 1);

        return vol;
    };

    /** Function that adds a disc portal.
     *
     * @param vol volume the portal should be added to
     * @param min_r lower radius of disc
     * @param max_r upper radius of disc
     * @param half_z z half length of the detector volume
     */
    auto add_disc_portal = [&](volume_type& vol, scalar min_r, scalar max_r,
                               scalar half_z) {
        m_id = {0, discs.size()};
        surfaces.emplace_back(transforms.size(), m_id, vol.index(),
                              inv_sf_finder);
        discs.emplace_back(min_r, max_r);
        // Position disc
        vector3 translation{0., 0., half_z};
        transforms.emplace_back(translation);

        vol.set_range({surfaces.size() - 1, surfaces.size()});
    };

    /** Function that adds a cylinder portal.
     *
     * @param vol volume the portal should be added to
     * @param r raius of the cylinder
     * @param half_z z half length of the cylinder
     */
    auto add_cylinder_portal = [&](volume_type& vol, scalar r, scalar half_z) {
        m_id = {1, cylinders.size()};
        surfaces.emplace_back(transforms.size(), m_id, vol.index(),
                              inv_sf_finder);
        cylinders.emplace_back(r, -half_z, half_z);
        transforms.emplace_back();  // identity

        vol.set_range({surfaces.size() - 1, surfaces.size()});
    };

    /** Function that updates surface and volume links when added to the global
     *  containers. The module surface edge point back into the volume.
     *
     * @param vol the volume the module surfaces belong to
     * @param modules the module surface container
     * @param trfs_offset offset into the global transform container
     * @param masks_offset offset into the global mask container
     */
    auto update_links = [&](volume_type& vol, surface_container& modules,
                            dindex trfs_offset, dindex masks_offset) {
        for (auto& sf : modules) {
            sf.transform() += trfs_offset;
            std::get<1>(sf.mask()) += masks_offset;
            sf.volume() = vol.index();
            sf.set_edge({vol.index(), inv_sf_finder});
        }

        vol.set_range({surfaces.size(), surfaces.size() + modules.size()});
    };

    //
    // beampipe
    //

    // build the beampipe volume
    volume_type& beampipe =
        add_cylinder_volume(0., beampipe_r, detector_half_z);

    // negative and positive, outer portal surface
    add_disc_portal(beampipe, 0., beampipe_r, -detector_half_z);
    add_disc_portal(beampipe, 0., beampipe_r, detector_half_z);
    add_cylinder_portal(beampipe, beampipe_r, detector_half_z);

    // Set disc portal edges
    dindex beampipe_portal_index = beampipe.range()[0];
    auto& beampipe_neg_portal = surfaces[beampipe_portal_index];
    beampipe_neg_portal.set_edge({leaving_world, inv_sf_finder});
    auto& beampipe_pos_portal = surfaces[++beampipe_portal_index];
    beampipe_pos_portal.set_edge({leaving_world, inv_sf_finder});

    //
    // first layer
    //

    // build the first layer volume
    volume_type& layer_1 =
        add_cylinder_volume(beampipe_r, first_layer_outer_r, detector_half_z);

    // negative and positive, inner and outer portal surface
    add_disc_portal(layer_1, beampipe_r, first_layer_outer_r, -detector_half_z);
    add_disc_portal(layer_1, beampipe_r, first_layer_outer_r, detector_half_z);
    add_cylinder_portal(layer_1, beampipe_r, detector_half_z);
    add_cylinder_portal(layer_1, first_layer_outer_r, detector_half_z);

    // Connect portals (beampipe, gap)
    const auto beampipe_idx = layer_1.index() - 1;
    auto& beampipe_outer_portal = surfaces[++beampipe_portal_index];
    beampipe_outer_portal.set_edge({layer_1.index(), inv_sf_finder});
    // discs
    dindex layer1_pt_index = layer_1.range()[0];
    auto& first_layer_neg_portal = surfaces[layer1_pt_index];
    first_layer_neg_portal.set_edge({leaving_world, inv_sf_finder});
    auto& first_layer_pos_portal = surfaces[++layer1_pt_index];
    first_layer_pos_portal.set_edge({leaving_world, inv_sf_finder});
    // cylinder
    auto& first_layer_inner_portal = surfaces[++layer1_pt_index];
    first_layer_inner_portal.set_edge({beampipe_idx, inv_sf_finder});

    // create module surfaces
    auto [l1_modules, l1_trfs, l1_masks] = create_modules<surface, rectangle>(
        8.4, 36., 0.145, 32., 2., 5., {16, 14});

    // update linking in volumes and surfaces
    update_links(layer_1, l1_modules, transforms.size(), rectangles.size());
    // Append to collections
    surfaces.insert(surfaces.end(), l1_modules.begin(), l1_modules.end());
    transforms.insert(transforms.end(), l1_trfs.begin(), l1_trfs.end());
    rectangles.insert(rectangles.end(), l1_masks.begin(), l1_masks.end());

    //
    // gap layer
    //

    // build gap volume
    volume_type& gap = add_cylinder_volume(
        first_layer_outer_r, second_layer_inner_r, detector_half_z);

    // negative and positive, inner and outer portal surface
    add_disc_portal(gap, first_layer_outer_r, second_layer_inner_r,
                    -detector_half_z);
    add_disc_portal(gap, first_layer_outer_r, second_layer_inner_r,
                    detector_half_z);
    add_cylinder_portal(gap, first_layer_outer_r, detector_half_z);
    add_cylinder_portal(gap, second_layer_inner_r, detector_half_z);

    // Connect portals (first layer, second layer)
    const auto layer_1_idx = gap.index() - 1;
    // Get the cylinder portals
    auto& first_layer_outer_portal = surfaces[++layer1_pt_index];
    first_layer_outer_portal.set_edge({gap.index(), inv_sf_finder});
    // discs
    dindex gap_pt_index = gap.range()[0];
    auto& gap_neg_portal = surfaces[gap_pt_index];
    gap_neg_portal.set_edge({leaving_world, inv_sf_finder});
    auto& gap_pos_portal = surfaces[++gap_pt_index];
    gap_pos_portal.set_edge({leaving_world, inv_sf_finder});
    // cylinder
    auto& gap_inner_portal = surfaces[++gap_pt_index];
    gap_inner_portal.set_edge({layer_1_idx, inv_sf_finder});

    //
    // second layer
    //

    // build the first layer volume
    volume_type& layer_2 = add_cylinder_volume(
        second_layer_inner_r, second_layer_outer_r, detector_half_z);

    // negative and positive, inner and outer portal surface
    add_disc_portal(layer_2, second_layer_inner_r, second_layer_outer_r,
                    -detector_half_z);
    add_disc_portal(layer_2, second_layer_inner_r, second_layer_outer_r,
                    detector_half_z);
    add_cylinder_portal(layer_2, second_layer_inner_r, detector_half_z);
    add_cylinder_portal(layer_2, second_layer_outer_r, detector_half_z);

    /// Connect portals (gap, world)
    const auto gap_idx = layer_2.index() - 1;
    // Get the cylinder portal
    auto& gap_outer_portal = surfaces[++gap_pt_index];
    gap_outer_portal.set_edge({layer_2.index(), inv_sf_finder});
    // discs
    dindex layer2_pt_index = layer_2.range()[0];
    auto& second_layer_neg_portal = surfaces[layer2_pt_index];
    second_layer_neg_portal.set_edge({leaving_world, inv_sf_finder});
    auto& second_layer_pos_portal = surfaces[++layer2_pt_index];
    second_layer_pos_portal.set_edge({leaving_world, inv_sf_finder});
    // cylinder
    auto& second_layer_inner_portal = surfaces[++layer2_pt_index];
    second_layer_inner_portal.set_edge({gap_idx, inv_sf_finder});
    auto& second_layer_outer_portal = surfaces[++layer2_pt_index];
    second_layer_outer_portal.set_edge({leaving_world, inv_sf_finder});

    // create module surfaces
    auto [l2_modules, l2_trfs, l2_masks] = create_modules<surface, rectangle>(
        8.4, 36., 0.145, 72., 2., 5., {32, 14});

    // update linking in volumes and surfaces
    update_links(layer_2, l2_modules, transforms.size(), rectangles.size());
    // Append to collections
    surfaces.insert(surfaces.end(), l2_modules.begin(), l2_modules.end());
    transforms.insert(transforms.end(), l2_trfs.begin(), l2_trfs.end());
    rectangles.insert(rectangles.end(), l2_masks.begin(), l2_masks.end());

    // Return all geometry containers
    return std::make_tuple<volume_container, surface_container,
                           transf_container, disc_container, cylinder_container,
                           rectangle_container>(
        std::move(volumes), std::move(surfaces), std::move(transforms),
        std::move(discs), std::move(cylinders), std::move(rectangles));
}
