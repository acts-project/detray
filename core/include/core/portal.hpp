/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "utils/containers.hpp"
#include "utils/enumerate.hpp"

#include <exception>
#include <memory>
#include <iostream>

namespace detray
{
    /** Portal class between volumes, it contains a surface/surface link,
     * and masks along and opposite the normal vector. As the number of
     * attached surfaces are usually very small, it only iterates through 
     * the corresponding masks and does not use more fancy logic.
     * 
     * This reduces the intersections to only one intersection even for
     * multiple associated volumes through one portal.
     * 
     * @tparam surface_link the surface or surface_link this mask represents
     * @param along_masks_t the type of mask container to be checked for an along side call
     * @param along_volume_offsets_t the type of volume offsets container to be checked for an along side call
     * @param opposite_masks_t the type of mask container to be checked for an opposite side call
     * @param opposite_volume_offsets_t the type of volume offsets container to be checked for an along side call
     */
    template <typename surface_link,
              typename along_masks_t,
              typename along_volume_offsets_t,
              typename opposite_masks_t = along_masks_t,
              typename opposite_volume_offsets_t = along_volume_offsets_t>
    struct portal
    {

        surface_link _surface;
        along_masks_t _along_masks;
        along_volume_offsets_t _along_offsets;
        opposite_masks_t _opposite_masks;
        opposite_volume_offsets_t _opposite_offsets;

        /** Contstuct a portal from arguments
         * 
         * @param surface the surface or surface link representation
         * @param along_masks the masks to be checked along the normal direction
         * @param along_offsets the offsets in the volume container for checks along the normal direction
         * @param opposite_masks the masks to be checked along the normal direction
         * @param opposite_offsets the offsets in the volume container for checks along the normal direction
         * 
         **/
        portal(surface_link &&surface, along_masks_t &&along_masks, along_volume_offsets_t &&along_offsets,
               opposite_masks_t &&opposite_masks, opposite_volume_offsets_t &&opposite_offsets)
            : _surface(std::move(surface)),
              _along_masks(std::move(along_masks)),
              _along_offsets(std::move(along_offsets)),
              _opposite_masks(std::move(opposite_masks)),
              _opposite_offsets(std::move(opposite_offsets))
        {
        }

        /** Access to the surface / surface link */
        const surface_link &surface() const { return _surface; }

        /** Get teh volume link along the normal direciton
         * 
         * @param p the point assumed to be on surface (i.e. local 2D/3D mask frame)
         * 
         **/
        template <typename point_type>
        const auto along(const point_type &p) const noexcept(false)
        {
            for (auto [i, mask] : enumerate(_along_masks))
            {
                if (mask(p) == e_inside)
                {
                    return i + _along_offsets[i];
                }
            }
            throw std::invalid_argument("portal: along volume requested, but point can not be associated.");
            return std::numeric_limits<decltype(_opposite_masks.size())>::max();
        }

        /** Get teh volume link opposite the normal direciton
         * 
         * @param p the point assumed to be on surface (i.e. local 2D/3D mask frame)
         * 
         **/
        template <typename point_type>
        const auto opposite(const point_type &p) const noexcept(false)
        {
            for (auto [i, mask] : enumerate(_opposite_masks))
            {
                if (mask(p) == e_inside)
                {
                    return i + _opposite_offsets[i];
                }
            }
            throw std::invalid_argument("portal: opposite volume requested, but point can not be associated.");
            return std::numeric_limits<decltype(_opposite_masks.size())>::max();
        }
    };

} // namespace detray
