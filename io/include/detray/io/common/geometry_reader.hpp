/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/io/common/detail/type_traits.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace std {

// Specialize std::hash so this pair can be used as key in e.g. a map
template <>
struct hash<std::pair<detray::surface_id, detray::io::detail::mask_shape>> {
    auto operator()(
        std::pair<detray::surface_id, detray::io::detail::mask_shape> p)
        const noexcept {
        // Calculate a hash from both numbers and avoid collisions
        return std::hash<std::size_t>()(static_cast<std::size_t>(p.first)) +
               1000u *
                   std::hash<std::size_t>()(static_cast<std::size_t>(p.second));
    }
};

}  // namespace std

namespace detray {

/// @brief Abstract base class for tracking geometry readers
template <class detector_t>
class geometry_reader : public reader_interface<detector_t> {

    using base_type = reader_interface<detector_t>;
    /// Can hold all types of surface fatory needed for the detector
    using sf_factory_ptr_t =
        std::shared_ptr<surface_factory_interface<detector_t>>;
    /// IO shape ids do not need to coincide with the detector mask ids,
    /// they are shared with ACTS
    using mask_shape = io::detail::mask_shape;
    /// Gets compile-time mask information
    template <mask_shape shape>
    using mask_info = detail::mask_info<shape, detector_t>;

    protected:
    /// Tag the reader as "geometry"
    inline static const std::string tag = "geometry";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Deserialize a detector @param det from its io payload @param det_data
    /// and add the volume names to @param name_map
    static void deserialize(detector_builder<typename detector_t::metadata,
                                             volume_builder>& det_builder,
                            typename detector_t::name_map& name_map,
                            const detector_payload& det_data) {

        // @todo Add volume grid

        // Deserialize the volumes one-by-one
        for (const auto& vol_data : det_data.volumes) {
            // Get a generic volume builder first and decorate it later
            auto vbuilder =
                det_builder.new_volume(static_cast<volume_id>(vol_data.type));

            // Set the volume name
            name_map[vbuilder->vol_index() + 1u] = vol_data.name;

            // @todo add the volume placement, once it can be checked for the
            // test detectors

            // Prepare the surface factories (one per shape and surface type)
            std::map<std::pair<surface_id, mask_shape>, sf_factory_ptr_t>
                sf_factories;

            // Add the surfaces to the factories
            for (const auto& sf_data : vol_data.surfaces) {

                const mask_payload& mask_data = sf_data.mask;

                // Check if a fitting factory already exists. If not, add it
                // dynamically
                const auto key = std::make_pair(sf_data.type, mask_data.shape);
                if (auto search = sf_factories.find(key);
                    search == sf_factories.end()) {
                    sf_factories[key] =
                        std::move(init_factory<mask_shape::n_shapes>(
                            mask_data.shape, sf_data.type));
                }

                // Add the data to the factory
                sf_factories.at(key)->push_back(deserialize(sf_data));
            }

            // Add the surfaces to the volume
            typename detector_t::geometry_context geo_ctx{};
            for (auto [key, sf_factory_ptr] : sf_factories) {
                // Sort the surfaces into the volume builder by type
                switch (sf_factory_ptr->surface_type()) {
                    case surface_id::e_portal:
                        vbuilder->add_portals(sf_factory_ptr, geo_ctx);
                        break;
                    case surface_id::e_sensitive:
                        vbuilder->add_sensitives(sf_factory_ptr, geo_ctx);
                        break;
                    case surface_id::e_passive:
                        vbuilder->add_passives(sf_factory_ptr, geo_ctx);
                        break;
                    case surface_id::e_unknown:
                        break;
                };
            }
        }

        det_builder.template set_volume_finder();
    }

    /// @returns a surface transform from its io payload @param trf_data
    static typename detector_t::transform3 deserialize(
        const transform_payload& trf_data) {
        using scalar_t = typename detector_t::scalar_type;
        using vector3_t = typename detector_t::vector3;

        vector3_t t{static_cast<scalar_t>(trf_data.tr[0]),
                    static_cast<scalar_t>(trf_data.tr[1]),
                    static_cast<scalar_t>(trf_data.tr[2])};
        vector3_t x{static_cast<scalar_t>(trf_data.rot[0]),
                    static_cast<scalar_t>(trf_data.rot[1]),
                    static_cast<scalar_t>(trf_data.rot[2])};
        vector3_t y{static_cast<scalar_t>(trf_data.rot[3]),
                    static_cast<scalar_t>(trf_data.rot[4]),
                    static_cast<scalar_t>(trf_data.rot[5])};
        vector3_t z{static_cast<scalar_t>(trf_data.rot[6]),
                    static_cast<scalar_t>(trf_data.rot[7]),
                    static_cast<scalar_t>(trf_data.rot[8])};

        return typename detector_t::transform3{t, x, y, z};
    }

    /// @returns surface data for a surface factory from a surface io payload
    /// @param trf_data
    static surface_data<detector_t> deserialize(
        const surface_payload& sf_data) {

        using nav_link_t = typename detector_t::surface_type::navigation_link;
        using scalar_t = typename detector_t::scalar_type;

        // Transcribe mask boundaries onto correct vector type
        std::vector<scalar_t> mask_boundaries;
        std::copy(sf_data.mask.boundaries.begin(),
                  sf_data.mask.boundaries.end(),
                  std::back_inserter(mask_boundaries));

        return {deserialize(sf_data.transform),
                static_cast<nav_link_t>(
                    base_type::deserialize(sf_data.mask.volume_link)),
                std::move(mask_boundaries)};
    }

    private:
    /// Determines the surface shape from the id @param shape_id in the payload
    /// and its type @param sf_type.
    ///
    /// @return the corresponding surface factory.
    template <mask_shape I>
    static sf_factory_ptr_t init_factory(
        const mask_shape shape_id, [[maybe_unused]] const surface_id sf_type) {
        // Shape index of surface data found
        if (shape_id == I) {
            // Get the corresponding mask shape type and its type id in the
            // detector
            using mask_info_t = mask_info<I>;
            using shape_t = typename mask_info_t::type;
            [[maybe_unused]] constexpr auto mask_id{mask_info_t::value};

            // Test wether this shape exists in detector
            if constexpr (not std::is_same_v<shape_t, void>) {
                // Get the correct factory for the type of surface
                switch (sf_type) {
                    case surface_id::e_portal:
                        using pt_factory_t =
                            surface_factory<detector_t, shape_t, mask_id,
                                            surface_id::e_portal>;
                        return std::make_shared<pt_factory_t>();
                    case surface_id::e_sensitive:
                        using sf_factory_t =
                            surface_factory<detector_t, shape_t, mask_id,
                                            surface_id::e_sensitive>;
                        return std::make_shared<sf_factory_t>();
                    case surface_id::e_passive:
                        using ps_factory_t =
                            surface_factory<detector_t, shape_t, mask_id,
                                            surface_id::e_passive>;
                        return std::make_shared<ps_factory_t>();
                    case surface_id::e_unknown:
                        throw std::runtime_error(
                            "Unknown surface type in geometry file");
                };
            }
        }
        // Test next shape id
        constexpr int current_id{static_cast<int>(I)};
        if constexpr (current_id > 0) {
            return init_factory<static_cast<mask_shape>(current_id - 1)>(
                shape_id, sf_type);
        }
        // Test some edge cases
        if (shape_id == mask_shape::unknown) {
            throw std::runtime_error("Unknown mask shape in geometry file!");
        } else {
            throw std::runtime_error(
                "Given mask type could not be matched: " +
                std::to_string(static_cast<std::int64_t>(shape_id)));

            // Cannot be reached
            return {nullptr};
        }
    }
};

}  // namespace detray
