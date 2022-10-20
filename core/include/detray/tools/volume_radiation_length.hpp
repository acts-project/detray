/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/utils/ranges.hpp"

// System include(s).
#include <climits>
#include <iostream>

namespace detray {

template <typename detector_t>
struct volume_radiation_length {

    using scalar_type = typename detector_t::scalar_type;

    struct material_property {
        scalar_type area = 0;
        scalar_type radiation_length = 0;
    };

    struct get_mask_area {
        using output_type = std::vector<scalar_type>;

        template <typename mask_group_t, typename index_t, typename surface_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const mask_group_t &mask_group, const index_t & /*index*/,
            const surface_t &surface) const {

            const auto &mask_range = surface.mask_range();

            output_type ret;

            for (const auto &mask :
                 detray::ranges::subrange(mask_group, mask_range)) {

                ret.push_back(mask.area());
            }
            return ret;
        }
    };

    struct get_material_property {

        using output_type = std::vector<material_property>;

        template <typename material_group_t, typename index_t,
                  typename surface_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const material_group_t &material_group, const index_t & /*index*/,
            const surface_t &surface) const {

            const auto &material_range = surface.material_range();

            output_type ret;

            for (const auto &mat :
                 detray::ranges::subrange(material_group, material_range)) {

                ret.push_back({mat.area(), mat.get_material().X0()});
            }

            return ret;
        }
    };

    volume_radiation_length(const detector_t &det) {

        // @todo: Consider the case where the volume is filled with gas
        const auto &mask_store = det.mask_store();
        const auto &mat_store = det.material_store();

        for (const auto &vol : det.volumes()) {
            const scalar_type vol_size = vol.volume_size();
            scalar_type denom = 0;

            for (const auto [obj_idx, obj] :
                 detray::views::enumerate(det.surfaces(), vol)) {

                const auto mask_areas =
                    mask_store.template call<get_mask_area>(obj.mask(), obj);

                const auto mat_properties =
                    mat_store.template call<get_material_property>(
                        obj.material(), obj);

                // Averaged radiation length (X_0) = V_0/sum(V_i/X_i) where
                // V_0 is the volume size and,
                // V_i and X_i is the size and radiation length of surface
                for (std::size_t i = 0; i < mask_areas.size(); i++) {
                    const scalar_type mat_size =
                        mask_areas[i] * mat_properties[i].area;

                    denom += mat_size / mat_properties[i].radiation_length;
                }
            }

            const scalar_type avg_rad_len = vol_size / denom;

            m_radiation_lengths.push_back(avg_rad_len);
        }
    }

    scalar_type operator[](std::size_t vol_index) const {
        return m_radiation_lengths[vol_index];
    }

    std::vector<scalar_type> operator()() const { return m_radiation_lengths; }

    private:
    std::vector<scalar_type> m_radiation_lengths = {};
};
}  // namespace detray