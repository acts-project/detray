/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <climits>
#include <map>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vector>

#include "detray/core/detector.hpp"
#include "detray/geometry/volume_connector.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/io/csv_io_types.hpp"
#include "detray/tools/bin_association.hpp"

namespace detray {

/// Function to read the detector from the CSV file
///
/// @tparam array_type is the type of the array used for the detector class
/// @tparam tuple_type is the type of the tuple used for the detector class
/// @tparam vector_type is the type of the tuple used for the detector class
///
/// All other container types are STL containers, as they do not leave the
/// function
///
/// @param detector_name is the name of the detector
/// @param surface_file_name is the name of the detector surface file
/// @param grid_file_name is the name of the surface grid file
/// @param layer_volume_file_name is the name of the file containing
/// layer/volume information
/// @param r_sync_tolerance is a tolerance to be for synching volumes in r
/// @param z_sync_tolerance is a toleranced to be used for synchinng volumes in
/// z
///
/// @return a detector object
template <typename detector_registry,
          template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple,
          template <typename...> class vector_type = dvector,
          template <typename...> class jagged_vector_type = djagged_vector,
          typename surface_source_link = dindex,
          typename bounds_source_link = dindex>
detector<detector_registry, array_type, tuple_type, vector_type,
         jagged_vector_type>
detector_from_csv(const std::string &detector_name,
                  const std::string &surface_file_name,
                  const std::string &layer_volume_file_name,
                  const std::string &grid_file_name,
                  const std::string &grid_entries_file_name,
                  std::map<dindex, std::string> &name_map,
                  vecmem::memory_resource &resource,
                  scalar /*r_sync_tolerance*/ = 0.,
                  scalar /*z_sync_tolerance*/ = 0.) {
    using alignable_store = static_transform_store<vector_type>;
    using detector_t = detector<detector_registry, array_type, tuple_type,
                                vector_type, jagged_vector_type>;

    name_map[0] = detector_name;
    detector_t d(resource);

    // Surface reading
    surface_reader s_reader(surface_file_name);
    csv_surface io_surface;

    // Layer volume reading
    layer_volume_reader lv_reader(layer_volume_file_name);
    csv_layer_volume io_layer_volume;

    // Surface grid reading
    surface_grid_reader sg_reader(grid_file_name);
    csv_surface_grid io_surface_grid;

    using volume_layer_index = std::pair<uint32_t, uint32_t>;
    std::map<volume_layer_index, typename detector_t::volume_type *> volumes;

    // Read in with a default context
    typename alignable_store::storage surface_transform_storage;
    typename alignable_store::context surface_default_context{};

    // Flushable containers
    typename detector_t::volume_type *c_volume = nullptr;
    typename detector_t::surface_filling_container c_surfaces;
    typename detector_t::mask_container c_masks(resource);
    typename detector_t::transform_container c_transforms;

    std::map<volume_layer_index, array_type<scalar, 6>> volume_bounds;

    // Remember the r/z attachements
    std::map<scalar, std::vector<dindex>> r_min_attachments;
    std::map<scalar, std::vector<dindex>> z_min_attachments;
    scalar r_max = 0;
    scalar z_max = -std::numeric_limits<scalar>::max();

    /** Helper method to attach volumes to bins
     *
     * @param attachments The attachnment map
     * @param value the (eventually new) value for insertion
     * @param volume_index the index of the attached volume
     */
    auto attach_volume = [](std::map<scalar, std::vector<dindex>> &attachments,
                            scalar value, dindex volume_index) -> void {
        if (attachments.find(value) == attachments.end()) {
            attachments[value] = {volume_index};
        } else {
            attachments[value].push_back(volume_index);
        }
    };

    // (A) Pre-read the bounds values
    std::map<scalar, scalar> z_min_layer_volumes;
    std::map<scalar, scalar> z_max_layer_volumes;
    std::map<scalar, scalar> r_min_layer_volumes;
    std::map<scalar, scalar> r_max_layer_volumes;

    while (lv_reader.read(io_layer_volume)) {
        volume_layer_index c_index = {io_layer_volume.volume_id,
                                      io_layer_volume.layer_id};
        array_type<scalar, 6> c_bounds = {
            io_layer_volume.min_v0, io_layer_volume.max_v0,
            io_layer_volume.min_v1, io_layer_volume.max_v1,
            io_layer_volume.min_v2, io_layer_volume.max_v2,
        };
        volume_bounds[c_index] = c_bounds;

        // Register low, high volume
        if (io_layer_volume.layer_id % 2 == 0) {
            z_min_layer_volumes[io_layer_volume.min_v1] =
                io_layer_volume.min_v1;
            z_max_layer_volumes[-1. * io_layer_volume.max_v1] =
                io_layer_volume.max_v1;
            r_min_layer_volumes[io_layer_volume.min_v0] =
                io_layer_volume.min_v0;
            r_max_layer_volumes[-1. * io_layer_volume.max_v0] =
                io_layer_volume.max_v0;
        }
    }
    /** Helper function to cluster boundaries
     *
     * @param boundaries is the unclustered boundaries map
     * @param tolerance is the tolerance parameter in which a cluster can lie
     * @param flip is a reverse flag for upper boundaries
     *
     */
    auto cluster_boundaries = [&](std::map<scalar, scalar> &boundaries,
                                  scalar tolerance, int flip = 1) -> void {
        scalar last_ = std::numeric_limits<scalar>::max();
        for (auto &[key, boundary] : boundaries) {
            // Do not adust the last one for max values
            if (flip < 0 and key == boundaries.begin()->first) {
                continue;
            }

            if (std::abs(last_ - flip * key) < tolerance) {
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
    auto find_and_replace = [](scalar &value,
                               const std::map<scalar, scalar> &value_map,
                               int flip = 1) -> void {
        // synchronize lower bound
        if (value_map.find(flip * value) != value_map.end()) {
            value = value_map.find(flip * value)->second;
        }
    };

    /** Helper function to return syncrhonized boundary objects
     *
     * @param bounds the unsynchronized bounds object
     * @param gap_volume an indicator if its a gap volume
     **/
    auto synchronize_bounds = [&](const array_type<scalar, 6> &bounds,
                                  bool gap_volume) -> array_type<scalar, 6> {
        scalar _r_min = bounds[0];
        scalar _r_max = bounds[1];
        scalar _z_min = bounds[2];
        scalar _z_max = bounds[3];
        scalar _phi_min = bounds[4];
        scalar _phi_max = bounds[5];

        if (not gap_volume) {
            find_and_replace(_r_min, r_min_layer_volumes);
            find_and_replace(_r_max, r_max_layer_volumes, -1);
            find_and_replace(_z_min, z_min_layer_volumes);
            find_and_replace(_z_max, z_max_layer_volumes, -1);
        } else {
            if (std::abs(_z_max + _z_min) < 0.1) {
                find_and_replace(_r_min, r_max_layer_volumes, -1);
                find_and_replace(_r_max, r_min_layer_volumes);
                find_and_replace(_z_min, z_min_layer_volumes);
                find_and_replace(_z_max, z_max_layer_volumes, -1);
            } else {
                find_and_replace(_r_min, r_min_layer_volumes);
                find_and_replace(_r_max, r_max_layer_volumes, -1);
                find_and_replace(_z_min, z_max_layer_volumes, -1);
                find_and_replace(_z_max, z_min_layer_volumes);
            }
        }

        return {_r_min, _r_max, _z_min, _z_max, _phi_min, _phi_max};
    };

    // Create the surface finders & reserve
    std::map<volume_layer_index, dindex> surface_finder_entries;

    using surfaces_r_axis = typename detector_t::surfaces_regular_axis;
    using surfaces_z_axis = typename detector_t::surfaces_regular_axis;
    using surfaces_phi_axis = typename detector_t::surfaces_circular_axis;

    using surfaces_r_phi_grid =
        typename detector_t::surfaces_regular_circular_grid;
    using surfaces_z_phi_grid =
        typename detector_t::surfaces_regular_circular_grid;
    using surfaces_finder = typename detector_t::surfaces_finder_type;

    surfaces_finder &detector_surfaces_finders = d.get_surfaces_finder();

    // (B) Pre-read the grids & create local object finders
    int sg_counts = 0;

    while (sg_reader.read(io_surface_grid)) {
        assert(sg_counts <
               static_cast<int>(detector_t::surfaces_finder_type::N_GRIDS));

        volume_layer_index c_index = {io_surface_grid.volume_id,
                                      io_surface_grid.layer_id};

        surface_finder_entries[c_index] = detector_surfaces_finders.size();

        bool is_disk = (io_surface_grid.type_loc0 == 3);

        // Prepare z axis parameters
        scalar _z_min = is_disk ? std::numeric_limits<scalar>::min()
                                : io_surface_grid.min_loc1;
        scalar _z_max = is_disk ? std::numeric_limits<scalar>::max()
                                : io_surface_grid.max_loc1;
        dindex _z_bins =
            is_disk ? 1u : static_cast<dindex>(io_surface_grid.nbins_loc1);

        // Prepare r axis parameters
        scalar _r_min = is_disk ? io_surface_grid.min_loc0 : 0.;
        scalar _r_max = is_disk ? io_surface_grid.max_loc0
                                : std::numeric_limits<scalar>::max();
        dindex _r_bins =
            is_disk ? static_cast<dindex>(io_surface_grid.nbins_loc0) : 1u;

        // Prepare phi axis parameters
        scalar _phi_min =
            is_disk ? io_surface_grid.min_loc1 : io_surface_grid.min_loc0;
        scalar _phi_max =
            is_disk ? io_surface_grid.max_loc1 : io_surface_grid.max_loc0;
        dindex _phi_bins =
            is_disk ? static_cast<dindex>(io_surface_grid.nbins_loc1)
                    : static_cast<dindex>(io_surface_grid.nbins_loc0);

        surfaces_z_axis z_axis{_z_bins, _z_min, _z_max, resource};
        surfaces_r_axis r_axis{_r_bins, _r_min, _r_max, resource};
        surfaces_phi_axis phi_axis{_phi_bins, _phi_min, _phi_max, resource};

        // negative / positive / inner / outer
        surfaces_r_phi_grid rphi_grid_n(r_axis, phi_axis, resource);
        surfaces_r_phi_grid rphi_grid_p(r_axis, phi_axis, resource);
        surfaces_z_phi_grid zphi_grid_i{z_axis, phi_axis, resource};
        surfaces_z_phi_grid zphi_grid_o{z_axis, phi_axis, resource};

        detector_surfaces_finders[sg_counts++] = std::move(rphi_grid_n);
        detector_surfaces_finders[sg_counts++] = std::move(rphi_grid_p);
        detector_surfaces_finders[sg_counts++] = std::move(zphi_grid_i);
        detector_surfaces_finders[sg_counts++] = std::move(zphi_grid_o);
    }

    typename detector_t::mask_defs::link_type mask_index = {dindex_invalid,
                                                            dindex_invalid};
    constexpr auto cylinder_id = detector_t::mask_defs::e_cylinder3;
    constexpr auto rectangle_id = detector_t::mask_defs::e_rectangle2;
    constexpr auto trapezoid_id = detector_t::mask_defs::e_trapezoid2;
    constexpr auto annulus_id = detector_t::mask_defs::e_annulus2;

    // (C) Read the surfaces and fill it
    while (s_reader.read(io_surface)) {
        volume_layer_index c_index = {io_surface.volume_id,
                                      io_surface.layer_id};
        auto c_volume_itr = volumes.find(c_index);
        if (c_volume_itr == volumes.end()) {
            // Flush the former information / c_volume still points to the prior
            // volume
            if (c_volume != nullptr) {
                d.add_objects(surface_default_context, *c_volume, c_surfaces,
                              c_masks, c_transforms);

                c_surfaces = typename detector_t::surface_filling_container();
                c_masks = typename detector_t::mask_container(resource);
                c_transforms = typename detector_t::transform_container();
            }

            // Find and fill the bounds
            auto new_bounds = volume_bounds.find(c_index);
            if (new_bounds == volume_bounds.end()) {
                // Bounds not found, do not build the volume
                continue;
            }

            const auto &unsynchronized_volume_bounds = new_bounds->second;
            // Check if you need to synchronize
            bool is_gap = (io_surface.layer_id % 2 != 0);
            auto _volume_bounds =
                synchronize_bounds(unsynchronized_volume_bounds, is_gap);

            // Check if this volume has a surface finder entry associated
            dindex surfaces_finder_entry = dindex_invalid;
            auto surface_finder_itr = surface_finder_entries.find(c_index);
            if (surface_finder_itr != surface_finder_entries.end()) {
                surfaces_finder_entry = surface_finder_itr->second;
            }

            std::string volume_name = detector_name;
            volume_name +=
                std::string("_vol_") + std::to_string(io_surface.volume_id);
            volume_name +=
                std::string("_lay_") + std::to_string(io_surface.layer_id);

            auto &new_volume =
                d.new_volume(_volume_bounds, surfaces_finder_entry);

            name_map[new_volume.index() + 1] = volume_name;

            // RZ attachment storage
            attach_volume(r_min_attachments, _volume_bounds[0],
                          new_volume.index());
            attach_volume(z_min_attachments, _volume_bounds[2],
                          new_volume.index());

            r_max = std::max(r_max, _volume_bounds[1]);
            z_max = std::max(z_max, _volume_bounds[3]);

            c_volume = &new_volume;
            // Insert to volume map
            volumes[c_index] = c_volume;
        } else {
            c_volume = c_volume_itr->second;
        }

        // Do not fill navigation layers
        if (io_surface.layer_id % 2 == 0) {

            // Read the transform
            vector3 t{io_surface.cx, io_surface.cy, io_surface.cz};
            vector3 x{io_surface.rot_xu, io_surface.rot_yu, io_surface.rot_zu};
            vector3 z{io_surface.rot_xw, io_surface.rot_yw, io_surface.rot_zw};

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
            if (bounds_type == 1) {
                // Cylinder Bounds

                // Add a new cylinder mask
                dindex cylinder_index = c_masks.template size<cylinder_id>();
                c_masks.template add_mask<cylinder_id>(
                    io_surface.bound_param0,
                    io_surface.cz - io_surface.bound_param1,
                    io_surface.cz + io_surface.bound_param1);
                // The read is valid: set the index
                mask_index = {cylinder_id, cylinder_index};

                // Build the cylinder transform
                auto &cylinder_transforms = std::get<cylinder_id>(c_transforms);
                cylinder_transforms.emplace_back(surface_default_context, t, z,
                                                 x);

                // Save the corresponding surface
                auto &cylinder_surfaces = c_surfaces[cylinder_id];
                cylinder_surfaces.emplace_back(
                    cylinder_transforms.size(surface_default_context) - 1,
                    mask_index, c_volume->index(), io_surface.geometry_id);
            } else if (bounds_type == 3) {
                // Disc bounds
            } else if (bounds_type == 6) {
                // Rectangle bounds

                // Add a new rectangle mask
                dindex rectangle_index = c_masks.template size<rectangle_id>();
                scalar half_x =
                    0.5 * (io_surface.bound_param2 - io_surface.bound_param0);
                scalar half_y =
                    0.5 * (io_surface.bound_param3 - io_surface.bound_param1);
                c_masks.template add_mask<rectangle_id>(half_x, half_y);
                // The read is valid: set the index
                mask_index = {rectangle_id, rectangle_index};

                // Build the rectangle transform
                auto &rectangle_transforms =
                    std::get<rectangle_id>(c_transforms);
                rectangle_transforms.emplace_back(surface_default_context, t, z,
                                                  x);

                // Save the corresponding surface
                auto &rectangle_surfaces = c_surfaces[rectangle_id];
                rectangle_surfaces.emplace_back(
                    rectangle_transforms.size(surface_default_context) - 1,
                    mask_index, c_volume->index(), io_surface.geometry_id);
            } else if (bounds_type == 7) {
                // Trapezoid bounds

                // Add a new trapezoid mask
                dindex trapezoid_index = c_masks.template size<trapezoid_id>();
                c_masks.template add_mask<trapezoid_id>(
                    io_surface.bound_param0, io_surface.bound_param1,
                    io_surface.bound_param2);

                // The read is valid: set the index
                mask_index = {trapezoid_id, trapezoid_index};

                // Build the trapezoid transform
                auto &trapezoid_transforms =
                    std::get<trapezoid_id>(c_transforms);
                trapezoid_transforms.emplace_back(surface_default_context, t, z,
                                                  x);

                // Save the corresponding surface
                auto &trapezoid_surfaces = c_surfaces[trapezoid_id];
                trapezoid_surfaces.push_back(
                    {trapezoid_transforms.size(surface_default_context) - 1,
                     mask_index, c_volume->index(), io_surface.geometry_id});
            } else if (bounds_type == 11) {
                // Annulus bounds

                // Add a new annulus mask
                dindex annulus_index = c_masks.template size<annulus_id>();
                c_masks.template add_mask<annulus_id>(
                    io_surface.bound_param0, io_surface.bound_param1,
                    io_surface.bound_param2, io_surface.bound_param3,
                    io_surface.bound_param4, io_surface.bound_param5,
                    io_surface.bound_param6);

                // The read is valid: set the index
                mask_index = {annulus_id, annulus_index};

                // Build the annulus transform
                auto &annulus_transforms = std::get<annulus_id>(c_transforms);
                annulus_transforms.emplace_back(surface_default_context, t, z,
                                                x);

                // Save the corresponding surface
                auto &annulus_surfaces = c_surfaces[annulus_id];
                annulus_surfaces.emplace_back(
                    annulus_transforms.size(surface_default_context) - 1,
                    mask_index, c_volume->index(), io_surface.geometry_id);
            }
        }  // end of exclusion for navigation layers
    }

    /** Helper method to sort and remove duplicates
     *
     * @param att attribute vector for sorting and duplicate removal
     *
     * @return the key values
     */
    auto sort_and_remove_duplicates =
        [](std::map<scalar, std::vector<dindex>> &att) -> dvector<scalar> {
        dvector<scalar> keys;
        keys.reserve(att.size());
        for (auto [key, value] : att) {
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
    axis::irregular<darray, dvector> raxis{{rs}, resource};
    axis::irregular<darray, dvector> zaxis{{zs}, resource};

    typename detector_t::volume_grid v_grid(std::move(raxis), std::move(zaxis),
                                            resource);

    // A step into the volume (stepsilon), can be read in from the smallest
    // difference
    scalar stepsilon = 1.;

    // Run the bin association and write out
    surface_grid_entries_writer sge_writer("grid-entries.csv");
    bool write_grid_entries =
        (grid_entries_file_name.find("write") != std::string::npos);
    bool read_grid_entries =
        not grid_entries_file_name.empty() and not write_grid_entries and
        not(grid_entries_file_name.find("none") != std::string::npos);

    // Loop over the volumes
    // - fill the volume grid
    // - run the bin association
    for (auto [iv, v] : enumerate(d.volumes())) {
        // Get the volume bounds for filling
        const auto &v_bounds = v.bounds();

        dindex irl = v_grid.axis_p0().bin(v_bounds[0] + stepsilon);
        dindex irh = v_grid.axis_p0().bin(v_bounds[1] - stepsilon);
        dindex izl = v_grid.axis_p1().bin(v_bounds[2] + stepsilon);
        dindex izh = v_grid.axis_p1().bin(v_bounds[3] - stepsilon);
        dindex volume_index = v.index();

        /*auto r_low = v_grid.axis_p0().borders(irl)[0];
        auto r_high = v_grid.axis_p0().borders(irh)[1];
        auto z_low = v_grid.axis_p1().borders(izl)[0];
        auto z_high = v_grid.axis_p1().borders(izh)[1];*/

        bool is_cylinder = std::abs(v_bounds[1] - v_bounds[0]) <
                           std::abs(v_bounds[3] - v_bounds[2]);

        for (dindex ir = irl; ir <= irh; ++ir) {
            for (dindex iz = izl; iz <= izh; ++iz) {
                v_grid.populate(ir, iz, std::move(volume_index));
            }
        }

        dindex sfi = v.surfaces_finder_entry();
        if (sfi != dindex_invalid and write_grid_entries) {
            auto &grid = is_cylinder ? detector_surfaces_finders[sfi + 2]
                                     : detector_surfaces_finders[sfi];
            bin_association(surface_default_context, d, v, grid, {0.1, 0.1},
                            false);

            csv_surface_grid_entry csv_ge;
            csv_ge.detray_volume_id = static_cast<int>(iv);
            size_t nbins0 = grid.axis_p0().bins();
            size_t nbins1 = grid.axis_p1().bins();
            for (size_t b0 = 0; b0 < nbins0; ++b0) {
                for (size_t b1 = 0; b1 < nbins1; ++b1) {
                    csv_ge.detray_bin0 = b0;
                    csv_ge.detray_bin1 = b1;
                    for (auto e : grid.bin(b0, b1)) {
                        csv_ge.detray_entry = e;
                        sge_writer.append(csv_ge);
                    }
                }
            }
        }
    }

    // Fast option, read the grid entries back in
    if (read_grid_entries) {

        surface_grid_entries_reader sge_reader(grid_entries_file_name);
        csv_surface_grid_entry surface_grid_entry;
        while (sge_reader.read(surface_grid_entry)) {
            // Get the volume bounds for fillind
            const auto &v =
                d.volume_by_index(surface_grid_entry.detray_volume_id);
            const auto &v_bounds = v.bounds();
            dindex sfi = v.surfaces_finder_entry();
            if (sfi != dindex_invalid) {
                bool is_cylinder = std::abs(v_bounds[1] - v_bounds[0]) <
                                   std::abs(v_bounds[3] - v_bounds[2]);
                auto &grid = is_cylinder ? detector_surfaces_finders[sfi + 2]
                                         : detector_surfaces_finders[sfi];
                // Fill the entry
                grid.populate(
                    static_cast<dindex>(surface_grid_entry.detray_bin0),
                    static_cast<dindex>(surface_grid_entry.detray_bin1),
                    static_cast<dindex>(surface_grid_entry.detray_entry));
            }
        }
    }

    // Connect the cylindrical volumes
    connect_cylindrical_volumes<detector_t, array_type, tuple_type,
                                vector_type>(d, v_grid);

    // Add the volume grid to the detector
    d.add_volume_grid(std::move(v_grid));

    return d;
}

}  // namespace detray
