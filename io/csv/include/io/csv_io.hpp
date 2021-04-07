/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "definitions/algebra.hpp"
#include "definitions/primitives.hpp"
#include "geometry/detector.hpp"
#include "geometry/grid.hpp"

#include <dfe/dfe_namedtuple.hpp>
#include <dfe/dfe_io_dsv.hpp>

#include <fstream>
#include <climits>
#include <map>

namespace detray
{

  struct csv_surface
  {
    /// Surface identifier. Not available in the TrackML datasets.
    uint64_t geometry_id;
    /// Partially decoded surface identifier components.
    uint32_t volume_id, boundary_id, layer_id, module_id;
    /// Center position components in mm.
    scalar cx, cy, cz;
    /// Rotation matrix components.
    scalar rot_xu, rot_xv, rot_xw;
    scalar rot_yu, rot_yv, rot_yw;
    scalar rot_zu, rot_zv, rot_zw;
    /// The type of the surface bpounds object, determines the parameters filled
    int bounds_type;
    scalar bound_param0 = -1.f;
    scalar bound_param1 = -1.f;
    scalar bound_param2 = -1.f;
    scalar bound_param3 = -1.f;
    scalar bound_param4 = -1.f;
    scalar bound_param5 = -1.f;
    scalar bound_param6 = -1.f;

    DFE_NAMEDTUPLE(csv_surface, geometry_id, volume_id, boundary_id, layer_id,
                   module_id, cx, cy, cz, rot_xu, rot_xv, rot_xw, rot_yu, rot_yv,
                   rot_yw, rot_zu, rot_zv, rot_zw, bounds_type, bound_param0,
                   bound_param1, bound_param2, bound_param3, bound_param4,
                   bound_param5, bound_param6);
  };

  using surface_reader = dfe::NamedTupleCsvReader<csv_surface>;

  struct csv_surface_grid
  {
    /// Surface identifier. Not available in the TrackML datasets.
    uint64_t geometry_id;
    /// Partially decoded surface identifier components.
    uint32_t volume_id, layer_id, surface_id;
    /// The number of bins in loc 0 / 1
    int type_loc0 = -1;
    int nbins_loc0 = -1;
    scalar min_loc0, max_loc0;
    int type_loc1 = -1;
    int nbins_loc1 = -1;
    scalar min_loc1, max_loc1;

    DFE_NAMEDTUPLE(csv_surface_grid, geometry_id, volume_id, layer_id, surface_id, type_loc1, nbins_loc0, min_loc0,
                   max_loc0, type_loc1, nbins_loc1, min_loc1, max_loc1);
  };

  using surface_grid_reader = dfe::NamedTupleCsvReader<csv_surface_grid>;


  /// Function to read the detector from the CSV file
  ///
  /// @param file_name the name of the file
  ///
  /// @return a drawable detector object
  inline static detector_map
  read_detector_map(const std::string &file_name, bool sensitives = true)
  {

/*
    detector c_detector;

    surface_reader reader(file_name);
    csv_surface io_surface;
    while (reader.read(io_surface))
    {
      // Check if we already have a volume, otherwise create one
      auto c_volume_itr = c_detector.find(io_surface.volume_id);
      if (c_volume_itr == c_detector.end()){
        volume n_volume;
        n_volume.volume_id = io_surface.volume_id;
        c_detector.insert({n_volume.volume_id, n_volume});
        c_volume_itr = c_detector.find(io_surface.volume_id);
      }
      // Check if we already have a layer, otherwise create one
      auto c_layer_itr = c_volume_itr->second.layers.find(io_surface.layer_id);
      if (c_layer_itr == c_volume_itr->second.layers.end()){
        layer n_layer;
        n_layer.volume_id = io_surface.volume_id;
        n_layer.layer_id = io_surface.layer_id;
        c_volume_itr->second.layers.insert({n_layer.layer_id, n_layer});
        c_layer_itr = c_volume_itr->second.layers.find(n_layer.layer_id);
      } 

      scalar layer_r_min = std::numeric_limits<scalar>::max();
      scalar layer_z_min = std::numeric_limits<scalar>::max();
      scalar layer_r_max = 0.;
      scalar layer_z_max = 0.;
      if (io_surface.module_id != 0 and sensitives)
      {
        aplot::surface surface;
        surface.geometry_id = io_surface.geometry_id;
        surface.volume_id = io_surface.volume_id;
        surface.layer_id = io_surface.layer_id;
        surface.module_id = io_surface.module_id;
        // Transform
        vector3 t{io_surface.cx, io_surface.cy, io_surface.cz};
        vector3 x{io_surface.rot_xu, io_surface.rot_yu, io_surface.rot_zu};
        vector3 z{io_surface.rot_xw, io_surface.rot_yw, io_surface.rot_zw};
        surface.transform = transform3{t, z, x};
        // Layer center positions for cylinder/disc detection
        scalar r = getter::perp(t);
        layer_r_min = std::min(layer_r_min, r);
        layer_z_min = std::min(layer_z_min, z[2]);
        layer_r_max = std::max(layer_r_max, r);
        layer_z_max = std::max(layer_z_min, z[2]);
        // bounds
        surface.bounds_type = io_surface.bounds_type;
        std::vector<scalar> bounds;
        bounds.push_back(io_surface.bound_param0);
        bounds.push_back(io_surface.bound_param1);
        bounds.push_back(io_surface.bound_param2);
        bounds.push_back(io_surface.bound_param3);
        bounds.push_back(io_surface.bound_param4);
        bounds.push_back(io_surface.bound_param5);
        bounds.push_back(io_surface.bound_param6);
        surface.bounds = bounds;
        // record it
        c_layer_itr->second.surfaces.push_back(std::move(surface));
      }
      // Detect cylinder 
      if ((layer_z_max-layer_z_min) > (layer_r_max-layer_r_min)){
        c_layer_itr->second.type = 0;
      } else {
        c_layer_itr->second.type = 1;
      }

    }
    return c_detector;
    */
  }


  /// Function to read the surface grid map from the CSV file
  ///
  /// @param file_name the name of the file
  ///
  /// @return a drawable detector object
  inline static grid_map 
  read_grid_map(const std::string& file_name){

  }


} // namespace aplot
