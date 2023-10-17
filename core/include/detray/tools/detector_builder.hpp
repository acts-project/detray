/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/definitions/bfield_backends.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/tools/grid_factory.hpp"
#include "detray/tools/volume_builder.hpp"
#include "detray/tools/volume_builder_interface.hpp"
#include "detray/utils/type_traits.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <memory>
#include <vector>

namespace detray {

/// @brief Provides functionality to build a detray detector volume by volume
///
/// @tparam metadata the type definitions for the detector
/// @tparam bfield_bknd_t the type of magnetic field to be used
/// @tparam volume_builder_t the basic volume builder to be used for the
///                          geometry data
/// @tparam volume_data_t the data structure that holds the volume builders
template <typename metadata = default_metadata,
          typename bfield_bknd_t = bfield::const_bknd_t,
          template <typename> class volume_builder_t = volume_builder,
          template <typename...> class volume_data_t = std::vector>
class detector_builder {
    public:
    using detector_type =
        detector<metadata, covfie::field<bfield_bknd_t>, host_container_types>;

    /// Empty detector builder
    detector_builder() = default;

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

    /// Decorate a volume builder at position @param volume_idx with more
    /// functionality
    template <class builder_t>
    DETRAY_HOST auto decorate(dindex volume_idx)
        -> volume_builder_interface<detector_type>* {

        m_volumes[volume_idx] =
            std::make_unique<builder_t>(std::move(m_volumes[volume_idx]));

        return m_volumes[volume_idx].get();
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

        detector_type det{resource};

        for (auto& vol_builder : m_volumes) {
            vol_builder->build(det);
        }

        det.set_bfield(std::move(m_bfield));
        det.set_volume_finder(std::move(m_vol_finder));

        // TODO: Add sorting, data deduplication etc. here later...

        return det;
    }

    /// Assembles the detector, without a magnetic field
    DETRAY_HOST
    void set_bfield(typename detector_type::bfield_type&& field) {
        m_bfield = std::forward<typename detector_type::bfield_type>(field);
    }

    /// Put the volumes into a search data structure
    template <typename... Args>
    DETRAY_HOST void set_volume_finder(Args&&... /*args*/) {

        using vol_finder_t = typename detector_type::volume_finder;

        // Add dummy volume grid for now
        if constexpr (detail::is_grid_v<vol_finder_t>) {

            // TODO: Construct it correctly with the grid builder
            mask<cylinder3D> vgrid_dims{0u,      0.f,   -constant<scalar>::pi,
                                        -2000.f, 180.f, constant<scalar>::pi,
                                        2000.f};
            std::array<std::size_t, 3> n_vgrid_bins{1u, 1u, 1u};

            grid_factory_type<vol_finder_t> vgrid_factory{};
            m_vol_finder = vgrid_factory.template new_grid<
                n_axis::open<n_axis::label::e_r>,
                n_axis::circular<n_axis::label::e_phi>,
                n_axis::open<n_axis::label::e_z>, n_axis::irregular<>,
                n_axis::regular<>, n_axis::irregular<>>(vgrid_dims,
                                                        n_vgrid_bins);
        }
    }

    protected:
    /// Data structure that holds a volume builder for every detector volume
    volume_data_t<std::unique_ptr<volume_builder_interface<detector_type>>>
        m_volumes{};
    /// Data structure to find volumes
    typename detector_type::volume_finder m_vol_finder{};
    /// The detector bfield
    typename detector_type::bfield_type m_bfield;
};

}  // namespace detray
