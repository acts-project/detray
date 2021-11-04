/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {
/// Connector method for cylindrical volumes without phi separation
///
/// @tparam detector_type is the type of the detector for the volume access
///
/// @param d [in,out] the detector to which the portal surfaces are added
/// @param volume_grid [in] the indexed volume grid
///
template <typename detector_type,
          template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple,
          template <typename> class vector_type = dvector>
void connect_cylindrical_volumes(
    detector_type &d, const typename detector_type::volume_grid &volume_grid) {
    typename detector_type::context default_context = {};

    // The grid is populated, now create portal surfaces
    // Start from left bottom corner (0,0)
    vector_type<array_type<dindex, 2>> seeds = {{0, 0}};
    dmap<dindex, dindex> seed_map;

    // The axes are used quite a bit
    const auto &axis_r = volume_grid.axis_p0();
    const auto &axis_z = volume_grid.axis_p1();

    /*** Helper function to add a new seed
     * It checks for validity of the seed
     *
     * @param seed the (potential new seed)
     * @param volume_index the volume index for this seed
     *        to avoid duplicate entries
     *
     * @note seeds are only set in bottom left corners of blocks
     **/
    auto add_new_seed = [&](const array_type<dindex, 2> &seed,
                            dindex volume_index) -> void {
        if (volume_index == dindex_invalid) {
            return;
        }

        if (seed_map.find(volume_index) != seed_map.end()) {
            ++seed_map[volume_index];
            return;
        }
        seed_map[volume_index] = 1;

        if (seed[0] < axis_r.bins() and seed[1] < axis_z.bins()) {
            seeds.push_back(seed);
        }
    };

    // Growing loop over the seeds, process to the end
    for (dindex iseed = 0; iseed < seeds.size(); ++iseed) {
        // The running seed & the refence seed
        auto seed = seeds[iseed];
        const auto &ref = volume_grid.bin(seed[0], seed[1]);

        // Build and add the portal surfaces
        auto &volume = d.volume_by_index(ref);

        // Collect portals per seed
        vector_type<tuple_type<array_type<scalar, 2>, dindex>>
            left_portals_info;
        vector_type<tuple_type<array_type<scalar, 2>, dindex>>
            upper_portals_info;
        vector_type<tuple_type<array_type<scalar, 2>, dindex>>
            right_portals_info;
        vector_type<tuple_type<array_type<scalar, 2>, dindex>>
            lower_portals_info;

        /// Helper method for walking up along the bins
        ///
        /// @param start_bin is the start bin of the walk
        /// @param portals_info is the container to collect for portals
        /// @param peek is the peek direction in z
        /// @param add_seed is a boolean whether new seeds should be added
        ///
        /// @return the end position of the the walk (inside position)
        auto walk_up =
            [&](array_type<dindex, 2> start_bin,
                vector_type<tuple_type<array_type<scalar, 2>, dindex>>
                    &portals_info,
                int peek, bool add_seed = false) -> array_type<dindex, 2> {
            auto running_bin = start_bin;
            array_type<dindex, 2> last_added = {dindex_invalid, dindex_invalid};
            // Test entry
            auto test = volume_grid.bin(running_bin[0], running_bin[1]);
            // Low/high
            auto low = axis_r.borders(running_bin[0]);
            auto high = axis_r.borders(running_bin[0]);
            //  1 - Left walk up
            // First peek to the to the left for portal desitnations
            dindex last_portal_dest =
                (running_bin[1] > 0 and running_bin[1] < axis_r.bins())
                    ? volume_grid.bin(running_bin[0], running_bin[1] + peek)
                    : dindex_invalid;
            while ((ref == test) and ++running_bin[0] < axis_r.bins()) {
                high = axis_r.borders(running_bin[0]);
                test = (seed[0] + 1 < axis_r.bins())
                           ? volume_grid.bin(running_bin[0], running_bin[1])
                           : dindex_invalid;
                // Peek outside and see if the portal destination has changed
                dindex portal_dest =
                    (running_bin[1] > 0 and running_bin[1] < axis_r.bins())
                        ? volume_grid.bin(running_bin[0], running_bin[1] + peek)
                        : dindex_invalid;
                if (portal_dest != last_portal_dest) {
                    // Record the boundary
                    portals_info.push_back(
                        {{low[0], high[0]}, last_portal_dest});
                    last_added = {running_bin[0] - 1, running_bin[1]};
                    // low is the new high
                    low = high;
                    last_portal_dest = portal_dest;
                }
                high = axis_r.borders(running_bin[0]);
            }

            // First Potential new seed
            if (add_seed) {
                add_new_seed(running_bin, test);
            }
            // By this we undo the overstepping in the loop (either by grid
            // boundary or ref/test fail)
            high = axis_r.borders(--running_bin[0]);
            if (last_added != running_bin) {
                portals_info.push_back({{low[0], high[1]}, last_portal_dest});
            }

            // The new position after walking is returned
            return running_bin;
        };

        /// Helper method for walking up along the bins
        ///
        /// @param start_bin is the start bin of the walk
        /// @param portals_info is the container to collect for portals
        /// @param peek is the peek direction in z
        /// @param add_seed is a boolean whether new seeds should be added
        /// @param walk_only is a boolean whether to actually add boundaries or
        /// not
        ///
        /// @return the end position of the the walk (inside position)
        auto walk_right =
            [&](const array_type<dindex, 2> &start_bin,
                vector_type<tuple_type<array_type<scalar, 2>, dindex>>
                    &portals_info,
                int peek, bool add_seed = false,
                bool walk_only = false) -> array_type<dindex, 2> {
            auto running_bin = start_bin;
            array_type<dindex, 2> last_added = {dindex_invalid, dindex_invalid};

            // Test, low and high at seed position
            auto test = volume_grid.bin(running_bin[0], running_bin[1]);
            // Low/high
            auto low = axis_z.borders(running_bin[1]);
            auto high = axis_z.borders(running_bin[1]);

            dindex last_portal_dest =
                (running_bin[0] < axis_r.bins())
                    ? volume_grid.bin(running_bin[0] + peek, running_bin[1])
                    : dindex_invalid;

            // Seed setting loop as well
            while (ref == test and (++running_bin[1]) < axis_z.bins()) {
                high = axis_z.borders(running_bin[1]);
                test = volume_grid.bin(running_bin[0], running_bin[1]);
                // Peek outside and see if the portal destination has changed
                dindex portal_dest =
                    (running_bin[0] < axis_r.bins())
                        ? volume_grid.bin(running_bin[0] + peek, running_bin[1])
                        : dindex_invalid;
                if (portal_dest != last_portal_dest) {
                    // Record the boundary
                    if (not walk_only) {
                        portals_info.push_back(
                            {{low[0], high[0]}, last_portal_dest});
                        last_added = {running_bin[0], running_bin[1] - 1};
                    }
                    // low is the new high
                    low = high;
                    last_portal_dest = portal_dest;
                    // That's a new seed right here, except for last one
                    if (running_bin[1] < axis_z.bins() and add_seed) {
                        add_new_seed({running_bin[0] + peek, running_bin[1]},
                                     portal_dest);
                    }
                }
                high = axis_z.borders(running_bin[1]);
            }
            // By this we undo the overstepping (see above)
            high = axis_z.borders(--running_bin[1]);
            if (not walk_only and last_added != running_bin) {
                portals_info.push_back({{low[0], high[1]}, last_portal_dest});
            }
            // The new bin position after walk
            return running_bin;
        };

        // Walk up from the (initial) seed position
        auto up_left = walk_up(seed, left_portals_info, true, -1);
        // Walk to the right from the resulting upper position
        walk_right(up_left, upper_portals_info, true, 1);
        // Walk right from the (initial) seed position
        auto bottom_right =
            walk_right(seed, lower_portals_info, false, -1, seed[0] == 0);
        // Walk up from the bottom right corner
        walk_up(bottom_right, right_portals_info, false, 1);

        typename detector_type::geometry::portal_filling_container portals = {};
        typename detector_type::mask_container portal_masks;
        typename detector_type::transform_container portal_transforms;

        // The bounds can be used for the mask and transform information
        const auto &volume_bounds = volume.bounds();

        /** Helper lambda to build disc portals
         *
         * @param portals_info The volume_index
         * @param bound_index The access for the boundary parameter
         *
         **/
        auto add_disc_portals =
            [&](vector_type<tuple_type<array_type<scalar, 2>, dindex>>
                    &portals_info,
                dindex bound_index) -> void {
            using portal_t = typename detector_type::geometry::portal;
            // Fill in the left side portals
            if (not portals_info.empty()) {
                // The portal transfrom is given from the left
                __plugin::vector3 _translation{0., 0.,
                                               volume_bounds[bound_index]};

                // Get the mask context group and fill it
                constexpr auto disc_id = detector_type::mask_id::e_portal_ring2;
                auto &disc_portal_transforms =
                    std::get<disc_id>(portal_transforms);
                auto &disc_portals = std::get<disc_id>(portals);

                /*typename detector_type::portal::mask_links mask_index = {
                    disc_id, portal_masks.template size<disc_id>()};
                // Create a stub mask for every unique index
                for (auto &info_ : portals_info) {
                    // Add new mask to container
                    portal_masks.template add_mask<disc_id>(
                        std::get<0>(info_)[0], std::get<0>(info_)[1]);

                    std::get<1>(mask_index) =
                        portal_masks.template size<disc_id>();

                    // Update mask link
                    portal_masks.template group<disc_id>().back().links() = {
                        std::get<1>(info_), dindex_invalid};

                    // Create the portal
                    typename detector_type::portal _portal{
                        disc_portal_transforms.size(default_context),
                        mask_index, volume.index(), dindex_invalid};
                    // Save the data
                    disc_portals.push_back(std::move(_portal));
                }*/

                typename portal_t::mask_links mask_index = {
                    disc_id,
                    {portal_masks.template size<disc_id>(),
                     portal_masks.template size<disc_id>()}};
                // Create a stub mask for every unique index
                for (auto &info_ : portals_info) {
                    // Add new mask to container
                    portal_masks.template add_mask<disc_id>(
                        std::get<0>(info_)[0], std::get<0>(info_)[1]);

                    std::get<1>(mask_index)[1] =
                        portal_masks.template size<disc_id>();

                    // Update mask link
                    portal_masks.template group<disc_id>().back().links() = {
                        std::get<1>(info_), dindex_invalid};
                }
                // Create the portal
                portal_t _portal{disc_portal_transforms.size(default_context),
                                 mask_index, volume.index(), dindex_invalid};
                // Save the data
                disc_portals.push_back(std::move(_portal));
                disc_portal_transforms.emplace_back(default_context,
                                                    _translation);
            }
        };

        /** Helper lambda for cylindrical portals
         *
         * @param portals_info
         * @param bound_index
         **/
        auto add_cylinder_portal =
            [&](vector_type<tuple_type<array_type<scalar, 2>, dindex>>
                    &portals_info,
                dindex bound_index) -> void {
            using portal_t = typename detector_type::geometry::portal;
            // Fill in the upper side portals
            if (not portals_info.empty()) {
                // Get the mask context group and fill it
                constexpr auto cylinder_id =
                    detector_type::mask_id::e_portal_cylinder3;
                auto &cylinder_portal_transforms =
                    std::get<cylinder_id>(portal_transforms);
                auto &cylinder_portals = std::get<cylinder_id>(portals);

                /*typename detector_type::portal::mask_links mask_index = {
                    cylinder_id, portal_masks.template size<cylinder_id>()};
                for (auto &info_ : portals_info) {
                    // Add new mask to container
                    const auto cylinder_range = std::get<0>(info_);

                    portal_masks.template add_mask<cylinder_id>(
                        volume_bounds[bound_index], cylinder_range[0],
                        cylinder_range[1]);

                    std::get<1>(mask_index) =
                        portal_masks.template size<cylinder_id>();

                    // Update masks links
                    portal_masks.template group<cylinder_id>().back().links() =
                        {std::get<1>(info_), dindex_invalid};

                    // Create the portal
                    typename detector_type::portal _portal{
                        cylinder_portal_transforms.size(default_context),
                        mask_index, volume.index(), dindex_invalid};
                    cylinder_portals.push_back(std::move(_portal));
                }*/

                typename portal_t::mask_links mask_index = {
                    cylinder_id,
                    {portal_masks.template size<cylinder_id>(),
                     portal_masks.template size<cylinder_id>()}};
                for (auto &info_ : portals_info) {
                    // Add new mask to container
                    const auto cylinder_range = std::get<0>(info_);

                    portal_masks.template add_mask<cylinder_id>(
                        volume_bounds[bound_index], cylinder_range[0],
                        cylinder_range[1]);

                    std::get<1>(mask_index)[1] =
                        portal_masks.template size<cylinder_id>();

                    // Update masks links
                    portal_masks.template group<cylinder_id>().back().links() =
                        {std::get<1>(info_), dindex_invalid};
                }
                // Create the portal
                portal_t _portal{
                    cylinder_portal_transforms.size(default_context),
                    mask_index, volume.index(), dindex_invalid};
                cylinder_portals.push_back(std::move(_portal));
                // This will be concentric targetted at nominal center
                cylinder_portal_transforms.emplace_back(default_context);
            }
        };

        // Add portals to the volume
        add_disc_portals(left_portals_info, 2);
        add_cylinder_portal(upper_portals_info, 1);
        add_disc_portals(right_portals_info, 3);
        add_cylinder_portal(lower_portals_info, 0);

        // Add portals to detector
        d.add_objects(default_context, volume, portals, portal_masks,
                      portal_transforms);
    }
}
}  // namespace detray
