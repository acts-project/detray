/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <cmath>
#include <map>
#include <string>

#include "core/detector.hpp"
#include "geometry/surface_base.hpp"
#include "geometry/volume.hpp"
#include "io/csv_io.hpp"
#include "masks/masks.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

namespace {

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
auto read_from_csv(detector_input_files& files) {
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

    auto d = detray::detector_from_csv<>(files.det_name, surfaces, volumes,
                                         grids, grid_entries, name_map);

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
auto create_modules(scalar m_half_x = 8.4, scalar m_half_y = 36.,
                    scalar m_tilt_phi = 0.145, scalar layer_r = 32.,
                    scalar radial_stagger = 2., scalar l_overlap = 5.,
                    const std::pair<int, int> binning = {16, 14}) {

    // using transform3 = __plugin::transform3;
    // using point3 = __plugin::point3;
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

    // First value is zero for rectangle type
    mask_index m_id = {0, 0};

    for (auto& m_center : m_centers) {

        // Surfaces with the linking into the local containers
        m_id = {0, masks.size()};
        surfaces.emplace_back(transforms.size(), m_id, detray::dindex_invalid,
                              detray::dindex_invalid);

        // The rectangle bounds for this module
        masks.emplace_back(m_half_x, m_half_y);

        // Build the transform
        // The local phi
        scalar m_phi = algebra::getter::phi(m_center);
        // Local z axis is the normal vector
        point3 m_local_z{std::cos(m_phi + m_tilt_phi),
                         std::sin(m_phi + m_tilt_phi), 0.};
        // Local x axis the normal to local y,z
        point3 m_local_x{-std::sin(m_phi + m_tilt_phi),
                         std::cos(m_phi + m_tilt_phi), 0.};

        // Create the module transform
        transforms.emplace_back(m_center, m_local_z, m_local_x);
    }

    return std::make_tuple<surface_container, trfs_container, mask_container>(
        std::move(surfaces), std::move(transforms), std::move(masks));
}

/** Add a single barrel layer volume to an existing collection.
 *
 * @param min_r minimal radius of volume
 * @param max_r maximal radius of volume
 * @param half_z half length in z of volume
 *
 * @return a detray cylinder volume
 */
template <typename volume_container_t>
auto& add_cylinder_volume(volume_container_t& volumes, scalar min_r = 25.,
                          scalar max_r = 40., scalar half_z = 500.) {

    // The volume bounds
    detray::darray<scalar, 6> bounds = {min_r,  max_r, -half_z,
                                        half_z, -M_PI, M_PI};

    // Add the new volume to the collection
    auto& v = volumes.emplace_back(bounds);
    v.set_index(volumes.size() - 1);

    return v;
}

/** Builds a simple detray geometry of the innermost tml layers. It contains:
 *
 * - a beampipe (r = 27mm, half_z = 500mm)
 * - a layer (r_min = 27mm, r_max = 38mm) with n rectangular modules at
 *   r = 32mm
 *
 * @returns a tuple containing the geometry objects collections: [volumes,
 *          surfaces, transforms, cylinder masks (portals), rectangle masks
 *          (modules)]
 */
auto toy_geometry() {
    constexpr bool for_surface = true;
    constexpr bool for_portal = false;

    // Volume type
    using volume_type = detray::volume<detray::darray>;
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

    // The surface transforms
    // using transform3 = __plugin::transform3;
    // We have cylinder and rectangle surfaces
    using cylinder = detray::cylinder3<false, detray::cylinder_intersector,
                                       __plugin::cylindrical2, edge_links, 0>;
    using rectangle = detray::rectangle2<detray::planar_intersector,
                                         __plugin::cartesian2, edge_links, 1>;

    // The surface type, both for volume portals and contained detector
    // surfaces
    using surface = detray::surface_base<transform_link, mask_index,
                                         volume_index, source_link, edge_links>;

    // The geometry data
    using volume_container = detray::dvector<volume_type>;
    using surface_container = detray::dvector<surface>;
    using transf_container = detray::dvector<transform3>;
    using cylinder_container = detray::dvector<cylinder>;
    using rectangle_container = detray::dvector<rectangle>;

    /** Volumes */
    volume_container volumes = {};
    /** Surfaces, including portals */
    surface_container surfaces = {};
    /** Surface transforms */
    transf_container transforms = {};
    /** Cylinder masks for volume boundaries (portals) */
    cylinder_container cylinders = {};
    /** Rectangle masks for detector surfaces */
    rectangle_container rectangles = {};
    // mask index for surfaces
    mask_index m_id = {};

    // parameters
    scalar detector_half_z = 500.;
    scalar beampipe_r = 27.;
    scalar first_layer_outer_r = 38.;
    scalar second_layer_inner_r = 64.;
    scalar second_layer_outer_r = 80.;

    /** Function that adds a cylinder portal with an identity transform to a
     * volume.
     */
    auto add_portal = [&](volume_type& vol, scalar r, scalar half_z) {
        cylinders.emplace_back(r, -half_z, half_z);
        transforms.emplace_back();  // identity
        m_id = {cylinders.size(), 0};
        surfaces.emplace_back(transforms.size(), m_id, vol.index(), 0);
        vol.template set_range<for_portal>(
            {surfaces.size() - 1, surfaces.size()});
    };

    /** Function that updates surface links when added to the global containers.
     */
    auto update_links = [&](volume_type& vol, surface_container& modules,
                            dindex trfs_offset, dindex masks_offset) {
        for (auto& sf : modules) {
            sf.transform() += trfs_offset;
            std::get<1>(sf.mask()) += masks_offset;
            sf.volume() = vol.index();
        }

        vol.template set_range<for_surface>(
            {surfaces.size(), surfaces.size() + modules.size()});
    };

    //
    // beampipe
    //

    // build the beampipe volume
    volume_type& v =
        add_cylinder_volume(volumes, 0., beampipe_r, detector_half_z);

    // portal surface to first layer
    add_portal(v, beampipe_r, detector_half_z);

    // module surfaces
    v.template set_range<for_surface>({0, 0});  // no modules

    //
    // first layer
    //

    // build the first layer volume
    v = add_cylinder_volume(volumes, beampipe_r, first_layer_outer_r,
                            detector_half_z);

    // inner and outer portal surface
    add_portal(v, beampipe_r, detector_half_z);
    add_portal(v, first_layer_outer_r, detector_half_z);

    // Connect it to the beampipe volume
    auto& beampipe_portal = surfaces.front();
    beampipe_portal.set_edge({v.index(), 0});
    auto& first_layer_inner_portal = surfaces[1];
    first_layer_inner_portal.set_edge({v.index() - 1, 0});

    // create module surfaces
    auto [l1_mods, l1_trfs, l1_masks] = create_modules<surface, rectangle>();

    // update linking
    update_links(v, l1_mods, transforms.size(), rectangles.size());
    // Append to collections
    surfaces.insert(surfaces.end(), l1_mods.begin(), l1_mods.end());
    transforms.insert(transforms.end(), l1_trfs.begin(), l1_trfs.end());
    rectangles.insert(rectangles.end(), l1_masks.begin(), l1_masks.end());

    //
    // first gap
    //

    // build gap volume
    v = add_cylinder_volume(volumes, first_layer_outer_r, second_layer_outer_r,
                            detector_half_z);

    // inner and outer portal surface
    add_portal(v, first_layer_outer_r, detector_half_z);
    add_portal(v, second_layer_outer_r, detector_half_z);

    // Connect it to the first layer
    dindex current_portal = v.template range<for_portal>()[0];
    auto& first_layer_outer_portal = surfaces[current_portal];
    first_layer_outer_portal.set_edge({v.index(), 0});
    auto& first_gap_inner_portal = surfaces[++current_portal];
    first_gap_inner_portal.set_edge({v.index() - 1, 0});

    /*
      // First gap
      ActsScalar secondLayerInnerR = 64.;
      auto firstGapBounds = std::make_unique<CylinderVolumeBounds>(
          firstLayerOuterR, secondLayerInnerR, detectorHalfZ);
      auto firstGap = DetectorVolume::makeShared(
          Transform3::Identity(), std::move(firstGapBounds), "BarrelGap0");

      // Second layer
      ActsScalar secondLayerOuterR = 80.;
      auto secondLayer = createBarrelVolume(secondLayerInnerR,
      secondLayerOuterR, detectorHalfZ, "BarrelLayer1", 8.4, 36.,
                                            0.145, 72., 2., 5., {32, 14});

      // The volumes in R
      std::vector<std::shared_ptr<DetectorVolume>> barrelVolumes = {
          beamPipe, firstLayer, firstGap, secondLayer};

      // Return the container in R
      return CylindricalContainerHelper::containerInR(std::move(barrelVolumes),
                                                      "BarrelWithTwoLayers");*/
    return std::make_tuple<volume_container, surface_container,
                           transf_container, cylinder_container,
                           rectangle_container>(
        std::move(volumes), std::move(surfaces), std::move(transforms),
        std::move(cylinders), std::move(rectangles));
}

}  // namespace