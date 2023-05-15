/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <tuple>

namespace detray {

/// @brief Bind transform and mask data together for surface building.
///
/// Surface data is kept in separate detector containers and linked by indices,
/// so during geometry building, great care has to be taken to make sure that
/// all components of a surface get sorted and linked into the containers
/// together and associated with the correct surface.
template <typename detector_t>
class surface_data {
    public:
    using navigation_link = typename detector_t::surface_type::navigation_link;

    DETRAY_HOST
    surface_data(
        const typename detector_t::transform3 &trf, navigation_link volume_link,
        const std::vector<typename detector_t::scalar_type> &mask_boundaries)
        : m_transform{trf},
          m_volume_link{volume_link},
          m_boundaries{mask_boundaries} {}

    DETRAY_HOST
    auto get_data()
        -> std::tuple<typename detector_t::transform3 &, navigation_link &,
                      std::vector<typename detector_t::scalar_type> &> {
        return std::tie(m_transform, m_volume_link, m_boundaries);
    }

    private:
    /// The surface placement
    typename detector_t::transform3 m_transform;
    /// The index of the volume that this surface links to
    navigation_link m_volume_link;
    // simple tuple of all mask types in the detector. Only one entry is filled
    // with the mask that corresponds to this specific surface.
    std::vector<typename detector_t::scalar_type> m_boundaries;
};

/// @brief How to generate surfaces with their corresponding masks and
/// transforms.
///
/// Can be hard coded surface generation or json reader.
template <typename detector_t>
class surface_factory_interface {
    public:
    using navigation_link = typename detector_t::surface_type::navigation_link;

    surface_factory_interface() = default;
    virtual ~surface_factory_interface() = default;

    /// @returns the number of surfaces the factory will produce
    DETRAY_HOST
    virtual dindex size() const = 0;

    /// @returns a surface id that corresponds to the surface type
    /// i.e. portal/sensisitve/passive
    DETRAY_HOST
    virtual surface_id surface_type() const = 0;

    /// Uniform interface to surface data from different sources
    DETRAY_HOST
    virtual void push_back(surface_data<detector_t> &&) = 0;

    DETRAY_HOST
    virtual auto push_back(std::vector<surface_data<detector_t>> &&)
        -> void /*std::tuple<typename detector_t::transform3 &, navigation_link
                   &, std::vector<typename detector_t::scalar_type> &>*/
        = 0;

    DETRAY_HOST
    virtual auto operator()(
        typename detector_t::volume_type &volume,
        typename detector_t::surface_container_t &surfaces,
        typename detector_t::transform_container &transforms,
        typename detector_t::mask_container &masks,
        typename detector_t::geometry_context ctx = {}) const
        -> dindex_range = 0;
};

}  // namespace detray
