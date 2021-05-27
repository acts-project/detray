/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/detector.hpp"
#include "core/volume_connector.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"

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

  struct csv_layer_volume
  {
    /// Surface identifier. Not available in the TrackML datasets.
    uint64_t geometry_id;
    /// Partially decoded surface identifier components.
    uint32_t volume_id, layer_id;
    /// The type of the surface bpounds object, determines the parameters filled
    int volume_type;
    float min_v0 = -1.f;
    float max_v0 = -1.f;
    float min_v1 = -1.f;
    float max_v1 = -1.f;
    float min_v2 = -1.f;
    float max_v2 = -1.f;

    DFE_NAMEDTUPLE(csv_layer_volume, geometry_id, volume_id, layer_id, min_v0,
                   max_v0, min_v1, max_v1, min_v2, max_v2);
  };

  using layer_volume_reader = dfe::NamedTupleCsvReader<csv_layer_volume>;

  /// Function to read the detector from the CSV file
  ///
  /// @param detector_name is the name of the detector
  /// @param surface_file_name is the name of the detector surface file
  /// @param grid_file_name is the name of the surface grid file
  ///
  /// @return a drawable detector object
  template <typename alignable_store = static_transform_store,
            typename surface_source_link = dindex,
            typename bounds_source_link = dindex>
  detector<alignable_store,
           surface_source_link,
           bounds_source_link>
  detector_from_csv(const std::string &detector_name,
                    const std::string &surface_file_name,
                    const std::string &grid_file_name,
                    const std::string &layer_volume_file_name,
                    scalar r_sync_tolerance = 0.,
                    scalar z_sync_tolerance = 0.)
  {

    detector d(detector_name);

    // Surface reading
    surface_reader s_reader(surface_file_name);
    csv_surface io_surface;

    // Surface grid reading
    surface_grid_reader sg_reader(grid_file_name);
    csv_surface_grid io_surface_grid;

    // Layer volume reading
    layer_volume_reader lv_reader(layer_volume_file_name);
    csv_layer_volume io_layer_volume;

    using volume_layer_index = std::pair<uint32_t, uint32_t>;
    std::map<volume_layer_index, typename detector<alignable_store, surface_source_link, bounds_source_link>::volume *> volumes;

    // Read in with a default context
    typename alignable_store::storage surface_transform_storage;
    typename alignable_store::context surface_default_context;

    // Flushable containers
    typename detector<alignable_store, surface_source_link, bounds_source_link>::volume *c_volume = nullptr;
    typename detector<alignable_store, surface_source_link, bounds_source_link>::surface_container c_surfaces;
    typename detector<alignable_store, surface_source_link, bounds_source_link>::surface_mask_container c_masks;

    std::map<volume_layer_index, darray<scalar, 6>> volume_bounds;

    // Remember the r/z attachements
    dmap<scalar, std::vector<dindex>> r_min_attachments;
    dmap<scalar, std::vector<dindex>> z_min_attachments;
    scalar r_max = 0;
    scalar z_max = -std::numeric_limits<scalar>::max();

    /** Helper method to attach volumes to bins
     * 
     * @param attachments The attachnment map
     * @param value the (eventually new) value for insertion
     * @param volume_index the index of the attached volume
     */
    auto attach_volume = [](dmap<scalar, std::vector<dindex>> &attachments, scalar value, dindex volume_index) -> void
    {
      if (attachments.find(value) == attachments.end())
      {
        attachments[value] = {volume_index};
      }
      else
      {
        attachments[value].push_back(volume_index);
      }
    };

    // Pre-read the bounds values
    dmap<scalar, scalar> z_min_layer_volumes;
    dmap<scalar, scalar> z_max_layer_volumes;
    dmap<scalar, scalar> r_min_layer_volumes;
    dmap<scalar, scalar> r_max_layer_volumes;

    while (lv_reader.read(io_layer_volume))
    {
      volume_layer_index c_index = {io_layer_volume.volume_id, io_layer_volume.layer_id};
      darray<scalar, 6> c_bounds = {
          io_layer_volume.min_v0,
          io_layer_volume.max_v0,
          io_layer_volume.min_v1,
          io_layer_volume.max_v1,
          io_layer_volume.min_v2,
          io_layer_volume.max_v2,
      };
      volume_bounds[c_index] = c_bounds;

      // Register low, high volume
      if (io_layer_volume.layer_id % 2 == 0)
      {
        z_min_layer_volumes[io_layer_volume.min_v1] = io_layer_volume.min_v1;
        z_max_layer_volumes[-1. * io_layer_volume.max_v1] = io_layer_volume.max_v1;
        r_min_layer_volumes[io_layer_volume.min_v0] = io_layer_volume.min_v0;
        r_max_layer_volumes[-1. * io_layer_volume.max_v0] = io_layer_volume.max_v0;
      }
    }
    /** Helper function to cluster boundaries 
     * 
     * @param boundaries is the unclustered boundaries map
     * @param tolerance is the tolerance parameter in which a cluster can lie
     * @param flip is a reverse flag for upper boundaries
     * 
     */
    auto cluster_boundaries = [&](dmap<scalar, scalar> &boundaries, scalar tolerance, int flip = 1) -> void
    {

      scalar last_ = std::numeric_limits<scalar>::max();
      for (auto &[key, boundary] : boundaries)
      {
        // Do not adust the last one for max values
        if (flip < 0 and key == boundaries.begin()->first){
          continue;
        }

        if (std::abs(last_ - flip * key) < tolerance)
        {
          boundary = last_;
        }
        last_ = boundary;
      }
    };

    // Cluster boundary synchronization
    cluster_boundaries(z_min_layer_volumes, 5.);
    cluster_boundaries(z_max_layer_volumes, 5., -1);
    cluster_boundaries(r_min_layer_volumes, 5.);
    cluster_boundaries(r_max_layer_volumes, 5., -1);

    /** Helper functions to find and replace in case
     * 
     * @param value is the value in question
     * @param value_map is the map for the replacement value
     * 
     **/
    auto find_and_replace = [](scalar &value, const dmap<scalar, scalar> &value_map, int flip = 1) -> void
    {
      // synchronize lower bound
      if (value_map.find(flip * value) != value_map.end())
      {
        value = value_map.find(flip * value)->second;
      }
    };

    /** Helper function to return syncrhonized boundary objects
     * 
     * @param bounds the unsynchronized bounds object
     * @param gap_volume an indicator if its a gap volume 
     **/
    auto synchronize_bounds = [&](const darray<scalar, 6> &bounds, bool gap_volume) -> darray<scalar, 6>
    {
      scalar r_min = bounds[0];
      scalar r_max = bounds[1];
      scalar z_min = bounds[2];
      scalar z_max = bounds[3];
      scalar phi_min = bounds[4];
      scalar phi_max = bounds[5];

      if (not gap_volume)
      {
        find_and_replace(r_min, r_min_layer_volumes);
        find_and_replace(r_max, r_max_layer_volumes, -1);
        find_and_replace(z_min, z_min_layer_volumes);
        find_and_replace(z_max, z_max_layer_volumes, -1);
      }
      else
      {
        if (std::abs(z_max + z_min) < 0.1)
        {
          find_and_replace(r_min, r_max_layer_volumes, -1);
          find_and_replace(r_max, r_min_layer_volumes);
          find_and_replace(z_min, z_min_layer_volumes);
          find_and_replace(z_max, z_max_layer_volumes, -1);
        }
        else
        {
          find_and_replace(r_min, r_min_layer_volumes);
          find_and_replace(r_max, r_max_layer_volumes, -1);
          find_and_replace(z_min, z_max_layer_volumes, -1);
          find_and_replace(z_max, z_min_layer_volumes);
        }
      }

      return {r_min, r_max, z_min, z_max, phi_min, phi_max};
    };

    // Reading the surfaces
    while (s_reader.read(io_surface))
    {
      volume_layer_index c_index = {io_surface.volume_id, io_surface.layer_id};
      auto c_volume_itr = volumes.find(c_index);
      if (c_volume_itr == volumes.end())
      {
        // Flush the former information
        if (c_volume != nullptr and not surface_transform_storage.empty())
        {
          // Construction with move semantics
          c_volume->add_surface_transforms(surface_default_context, std::move(surface_transform_storage));
          c_volume->add_surface_components(std::move(c_surfaces), std::move(c_masks));
          // Get new clean containers
          surface_transform_storage = typename alignable_store::storage();
          c_surfaces =
              typename detector<alignable_store, surface_source_link, bounds_source_link>::surface_container();
          c_masks =
              typename detector<alignable_store, surface_source_link, bounds_source_link>::surface_mask_container();
        }

        // Create a new volume & assign
        std::string volume_name = detector_name;
        volume_name += std::string("_vol_") + std::to_string(io_surface.volume_id);
        volume_name += std::string("_lay_") + std::to_string(io_surface.layer_id);
        // Find and fill the bounds
        auto new_bounds = volume_bounds.find(c_index);
        if (new_bounds == volume_bounds.end())
        {
          // Bounds not found, do not build the volume
          continue;
        }

        const auto &unsynchronized_volume_bounds = new_bounds->second;
        // Check if you need to synchronize
        bool is_gap = (io_surface.layer_id % 2 != 0);
        auto volume_bounds = synchronize_bounds(unsynchronized_volume_bounds, is_gap);

        auto &new_volume = d.new_volume(volume_name, volume_bounds);

        // RZ attachment storage
        attach_volume(r_min_attachments, volume_bounds[0], new_volume.index());
        attach_volume(z_min_attachments, volume_bounds[2], new_volume.index());

        r_max = std::max(r_max, volume_bounds[1]);
        z_max = std::max(z_max, volume_bounds[3]);

        c_volume = &new_volume;
        // Insert to volume map
        volumes[c_index] = c_volume;
      }
      else
      {
        c_volume = c_volume_itr->second;
      }

      // Do not fill navigation layers
      if (io_surface.layer_id % 2 == 0)
      {

        // Read the transform
        vector3 t{io_surface.cx, io_surface.cy, io_surface.cz};
        vector3 x{io_surface.rot_xu, io_surface.rot_yu, io_surface.rot_zu};
        vector3 z{io_surface.rot_xw, io_surface.rot_yw, io_surface.rot_zw};
        dindex transform_index = surface_transform_storage.size();
        surface_transform_storage.push_back(transform3{t, z, x});

        // Translate the mask & add it to the mask container
        unsigned int bounds_type = io_surface.bounds_type;
        std::vector<scalar> bounds;
        bounds.push_back(io_surface.bound_param0);
        bounds.push_back(io_surface.bound_param1);
        bounds.push_back(io_surface.bound_param2);
        bounds.push_back(io_surface.bound_param3);
        bounds.push_back(io_surface.bound_param4);
        bounds.push_back(io_surface.bound_param5);
        bounds.push_back(io_surface.bound_param6);

        // Acts naming convention for bounds
        typename detector<alignable_store, surface_source_link, bounds_source_link>::surface_mask_index
            mask_index = {dindex_invalid, dindex_invalid};

        if (bounds_type == 1)
        {
          // Cylinder bounds
          constexpr auto cylinder_context = detector<alignable_store,
                                                     surface_source_link,
                                                     bounds_source_link>::surface_cylinder::mask_context;

          // Get the cylinder mask container
          auto &cylinder_masks = std::get<cylinder_context>(c_masks);
          dindex cylinder_index = cylinder_masks.size();
          cylinder_masks.push_back({io_surface.bound_param0, io_surface.cz - io_surface.bound_param1, io_surface.cz + io_surface.bound_param1});
          // The read is valid: set the index
          mask_index = {cylinder_context, cylinder_index};
        }
        else if (bounds_type == 3)
        {
          // Disc bounds
        }
        else if (bounds_type == 6)
        {
          // Rectangle bounds
          constexpr auto rectangle_context = detector<alignable_store,
                                                      surface_source_link,
                                                      bounds_source_link>::surface_rectangle::mask_context;

          // Get the rectangle mask container
          auto &rectangle_masks = std::get<rectangle_context>(c_masks);
          dindex rectangle_index = rectangle_masks.size();
          scalar half_x = 0.5 * (io_surface.bound_param2 - io_surface.bound_param0);
          scalar half_y = 0.5 * (io_surface.bound_param3 - io_surface.bound_param1);
          rectangle_masks.push_back({half_x, half_y});
          // The read is valid: set the index
          mask_index = {rectangle_context, rectangle_index};
        }
        else if (bounds_type == 7)
        {
          // Trapezoid bounds
          constexpr auto trapezoid_context = detector<alignable_store,
                                                      surface_source_link,
                                                      bounds_source_link>::surface_trapezoid::mask_context;
          // Get the trapezoid mask container
          auto &trapezoid_masks = std::get<trapezoid_context>(c_masks);
          dindex trapezoid_index = trapezoid_masks.size();
          trapezoid_masks.push_back({io_surface.bound_param0, io_surface.bound_param1, io_surface.bound_param2});
          // The read is valid: set the index
          mask_index = {trapezoid_context, trapezoid_index};
        }
        else if (bounds_type == 11)
        {
          // Annulus bounds
          constexpr auto annulus_context = detector<alignable_store,
                                                    surface_source_link,
                                                    bounds_source_link>::surface_annulus::mask_context;
          // Get the trapezoid mask container
          auto &annulus_masks = std::get<annulus_context>(c_masks);
          dindex annulus_index = annulus_masks.size();
          annulus_masks.push_back({io_surface.bound_param0,
                                   io_surface.bound_param1,
                                   io_surface.bound_param2,
                                   io_surface.bound_param3,
                                   io_surface.bound_param4,
                                   io_surface.bound_param5,
                                   io_surface.bound_param6});
          // The read is valid: set the index
          mask_index = {annulus_context, annulus_index};
        }

        // Fill the surface into the temporary container
        if (mask_index[0] != dindex_invalid)
        {
          c_surfaces.push_back({transform_index, mask_index, c_volume->index(), io_surface.geometry_id});
        }
      } // end of exclusion for navigation layers
    }

    /** Helper method to sort and remove duplicates
     * 
     * @param att attribute vector for sorting and duplicate removal
     * 
     * @return the key values
     */
    auto
        sort_and_remove_duplicates = [](dmap<scalar, std::vector<dindex>> &att) -> dvector<scalar>
    {
      dvector<scalar> keys;
      keys.reserve(att.size());
      for (auto [key, value] : att)
      {
        keys.push_back(key);
        std::sort(value.begin(), value.end());
        value.erase(std::unique(value.begin(), value.end()), value.end());
      }
      return keys;
    };

    // Drawing the lines for the grid search
    auto rs = sort_and_remove_duplicates(r_min_attachments);
    rs.push_back(r_max);
    auto zs = sort_and_remove_duplicates(z_min_attachments);
    zs.push_back(z_max);

    // Create axes and volume grid
    axis::irregular raxis{{rs}};
    axis::irregular zaxis{{zs}};

    typename detector<alignable_store, surface_source_link, bounds_source_link>::volume_grid
        v_grid(std::move(raxis), std::move(zaxis));

    // A step into the volume (stepsilon), can be read in from the smallest difference
    scalar stepsilon = 1.;

    // Loop over the volumes and fill the volume grid
    for (const auto &v : d.volumes())
    {
      // Get the volume bounds for fillind
      const auto &v_bounds = v.bounds();

      dindex irl = v_grid.axis_p0().bin(v_bounds[0] + stepsilon);
      dindex irh = v_grid.axis_p0().bin(v_bounds[1] - stepsilon);
      dindex izl = v_grid.axis_p1().bin(v_bounds[2] + stepsilon);
      dindex izh = v_grid.axis_p1().bin(v_bounds[3] - stepsilon);
      dindex volume_index = v.index();

      auto r_low = v_grid.axis_p0().borders(irl)[0];
      auto r_high = v_grid.axis_p0().borders(irh)[1];
      auto z_low = v_grid.axis_p1().borders(izl)[0];
      auto z_high = v_grid.axis_p1().borders(izh)[1];

      for (dindex ir = irl; ir <= irh; ++ir)
      {
        for (dindex iz = izl; iz <= izh; ++iz)
        {
          v_grid.populate(ir, iz, std::move(volume_index));
        }
      }
    }

    // Connect the cylindrical volumes
    connect_cylindrical_volumes(d, v_grid);

    // Add the volume grid to the detector
    d.add_volume_grid(std::move(v_grid));

    return d;
  }

} // namespace detray
