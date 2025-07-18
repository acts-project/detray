/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/builders/grid_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/utils/grid/detail/concepts.hpp"
#include "detray/utils/type_traits.hpp"

// Vecmem include(s)
#include <detray/utils/log.hpp>
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace detray {

/// @brief Provides functionality to build a detray detector volume by volume
///
/// @tparam metadata the type definitions for the detector
/// @tparam bfield_bknd_t the type of magnetic field to be used
/// @tparam volume_builder_t the basic volume builder to be used for the
///                          geometry data
/// @tparam volume_data_t the data structure that holds the volume builders
template <typename metadata,
          template <typename> class volume_builder_t = volume_builder,
          template <typename...> class volume_data_t = std::vector>
class detector_builder {
    public:
    using detector_type = detector<metadata, host_container_types>;
    using algebra_type = typename detector_type::algebra_type;
    using scalar_type = dscalar<algebra_type>;

    /// Set the name of the detector under construction to @param det_name
    DETRAY_HOST void set_name(std::string det_name) {
        m_detector_name = std::move(det_name);
    }

    /// @returns the name of the detector under construction
    DETRAY_HOST std::string_view name() const { return m_detector_name; }

    /// Add a new volume builder that will build a volume of the shape given by
    /// @param id
    template <typename... Args>
    DETRAY_HOST auto new_volume(const volume_id id, Args&&... args)
        -> volume_builder_interface<detector_type>* {

        m_volumes.push_back(std::make_unique<volume_builder_t<detector_type>>(
            id, static_cast<dindex>(m_volumes.size()),
            std::forward<Args>(args)...));

        return m_volumes.back().get();
    }

    /// @returns the number of volumes currently registered in the builder
    DETRAY_HOST auto n_volumes() const -> dindex {
        return static_cast<dindex>(m_volumes.size());
    }

    /// @returns 'true' if there is a volume builder registered for
    /// the volume with index @param volume_idx
    DETRAY_HOST bool has_volume(const std::size_t volume_idx) const {
        return volume_idx < m_volumes.size();
    }

    /// Decorate a volume builder at position @param volume_idx with more
    /// functionality
    template <class builder_t>
    DETRAY_HOST auto decorate(dindex volume_idx) -> builder_t* {
        assert(has_volume(volume_idx));

        m_volumes[volume_idx] =
            std::make_unique<builder_t>(std::move(m_volumes[volume_idx]));

        // Always works, we set it as this type in the line above
        return dynamic_cast<builder_t*>(m_volumes[volume_idx].get());
    }

    /// Decorate a volume builder @param v_builder with more functionality
    template <class builder_t>
    DETRAY_HOST auto decorate(const volume_builder_interface<detector_type>*
                                  v_builder) -> builder_t* {
        assert(v_builder != nullptr);

        return decorate<builder_t>(v_builder->vol_index());
    }

    /// Access a particular volume builder by volume index @param volume_idx
    DETRAY_HOST
    auto operator[](dindex volume_idx)
        -> volume_builder_interface<detector_type>* {
        return m_volumes[volume_idx].get();
    }

    /// Assembles the final detector from the volumes builders and allocates
    /// the detector containers with the memory resource @param resource
    DETRAY_HOST
    auto build(vecmem::memory_resource& resource) -> detector_type {

        DETRAY_DEBUG("detray: building detector " << name());

        detector_type det{resource};

        DETRAY_DEBUG("Have " << m_volumes.size()
                             << " configured volume builders");
        for (auto& vol_builder : m_volumes) {

            DETRAY_DEBUG("- builder: " << vol_builder->name());
            vol_builder->build(det);
        }

        DETRAY_DEBUG("Setting volume finder");
        det.set_volume_finder(std::move(m_vol_finder));

        // TODO: Add sorting, data deduplication etc. here later...

        DETRAY_DEBUG("detray: detector building complete");
        return det;
    }

    /// Assembles the final detector and fill an externally provided name map
    /// @param name_map
    DETRAY_HOST
    auto build(vecmem::memory_resource& resource,
               typename detector_type::name_map& name_map) -> detector_type {

        DETRAY_DEBUG("detray: filling names for detector " << name());

        assert(name_map.empty());

        // By convention the name of the detector is at position 0
        name_map.emplace(0u, m_detector_name);

        for (auto& vol_builder : m_volumes) {
            name_map.emplace(vol_builder->vol_index() + 1u,
                             vol_builder->name());
        }

        return build(resource);
    }

    /// Put the volumes into a search data structure
    template <typename... Args>
    DETRAY_HOST void set_volume_finder([[maybe_unused]] Args&&... args) {
        DETRAY_DEBUG("Setting volume finder for detector " << name());

        using vol_finder_t = typename detector_type::volume_finder;

        // Add dummy volume grid for now
        if constexpr (concepts::grid<vol_finder_t>) {

            // TODO: Construct it correctly with the grid builder
            mask<cylinder3D, algebra_type> vgrid_dims{
                0u,      0.f,   -constant<scalar_type>::pi,
                -2000.f, 180.f, constant<scalar_type>::pi,
                2000.f};
            darray<std::size_t, 3> n_vgrid_bins{1u, 1u, 1u};

            darray<std::vector<scalar_type>, 3UL> bin_edges{
                std::vector<scalar_type>{0.f, 180.f},
                std::vector<scalar_type>{-constant<scalar_type>::pi,
                                         constant<scalar_type>::pi},
                std::vector<scalar_type>{-2000.f, 2000.f}};

            grid_factory_type<vol_finder_t> vgrid_factory{};
            m_vol_finder = vgrid_factory.template new_grid<
                axis::open<axis::label::e_r>,
                axis::circular<axis::label::e_phi>,
                axis::open<axis::label::e_z>, axis::irregular<scalar_type>,
                axis::regular<scalar_type>, axis::irregular<scalar_type>>(
                vgrid_dims, n_vgrid_bins, {}, bin_edges);
        } else {
            m_vol_finder = vol_finder_t{args...};
        }
    }

    /// @returns access to the volume finder
    DETRAY_HOST typename detector_type::volume_finder& volume_finder() {
        return m_vol_finder;
    }

    private:
    /// Name of the new detector
    std::string m_detector_name{"detray_detector"};
    /// Data structure that holds a volume builder for every detector volume
    volume_data_t<std::unique_ptr<volume_builder_interface<detector_type>>>
        m_volumes{};
    /// Data structure to find volumes
    typename detector_type::volume_finder m_vol_finder{};
};

}  // namespace detray
