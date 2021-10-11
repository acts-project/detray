/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <iostream>
#include <map>
#include <string>

#include "core/detector.hpp"
#include "geometry/surface_base.hpp"
#include "geometry/volume.hpp"
#include "io/csv_io.hpp"
#include "masks/masks.hpp"

/// @note __plugin has to be defined with a preprocessor command
namespace detray_tests {

/** Keeps the relevant csv file names */
struct detector_input_files {
    std::string det_name, surface, layer_volume, surface_grid,
        surface_grid_entries;
};

/// open data detector
detector_input_files odd_files = {"odd", "odd.csv", "odd-layer-volumes.csv",
                                  "odd-surface-grids.csv", ""};

/// track ml detector
detector_input_files tml_files = {"tml", "tml.csv", "tml-layer-volumes.csv",
                                  "tml-surface-grids.csv", ""};

/** Read a detector from csv files */
auto read_from_csv(detector_input_files &files) {
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

// Build toy geometry
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
    using transform3 = __plugin::transform3;
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

    mask_index m_id = {};
    detray::darray<detray::scalar, 6> bounds = {};

    // beampipe
    detray::scalar detectorHalfZ = 500.;
    detray::scalar beamPipeR = 27.;

    cylinders.emplace_back(beamPipeR, detectorHalfZ, detectorHalfZ);
    transforms.emplace_back();  // identity
    m_id = {0, 0};
    surfaces.emplace_back(0, m_id, 0, 0);  // only one cylinder portal
    // build the volume
    bounds = {0, beamPipeR, -detectorHalfZ, detectorHalfZ, -M_PI, M_PI};
    volume_type &v = volumes.emplace_back(bounds);
    v.set_index(volumes.size() - 1);
    v.template set_range<for_portal>({0, 1});
    v.template set_range<for_surface>({0, 0});
    v.set_surfaces_finder(0);

    /*std::shared_ptr<DetectorVolume> createCentralDetector(
    ActsScalar detectorHalfZ = 500.) {
  // Create the volume bounds
  ActsScalar beamPipeR = 27.;
  // Beam pipe volume
  auto beamPipeBounds =
      std::make_unique<CylinderVolumeBounds>(0., beamPipeR, detectorHalfZ);
  auto beamPipe = DetectorVolume::makeShared(
      Transform3::Identity(), std::move(beamPipeBounds), "BeamPipe");
  // First layer
  ActsScalar firstLayerOuterR = 38.;
  auto firstLayer = createBarrelVolume(beamPipeR, firstLayerOuterR,
                                       detectorHalfZ, "BarrelLayer0");
  // First gap
  ActsScalar secondLayerInnerR = 64.;
  auto firstGapBounds = std::make_unique<CylinderVolumeBounds>(
      firstLayerOuterR, secondLayerInnerR, detectorHalfZ);
  auto firstGap = DetectorVolume::makeShared(
      Transform3::Identity(), std::move(firstGapBounds), "BarrelGap0");
  // Second layer
  ActsScalar secondLayerOuterR = 80.;
  auto secondLayer = createBarrelVolume(secondLayerInnerR, secondLayerOuterR,
                                        detectorHalfZ, "BarrelLayer1", 8.4, 36.,
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

}  // namespace detray_tests