/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <tuple>

namespace detray {

/// @brief How to generate surfaces with their corresponding masks and
/// transforms.
///
/// Can be hard coded surface generation or json reader.
template <typename detector_t>
class surface_factory_interface {
    public:
    surface_factory_interface() = default;
    virtual ~surface_factory_interface() = default;

    DETRAY_HOST
    virtual auto operator()(
        const typename detector_t::volume_type &volume,
        typename detector_t::surface_container &surfaces,
        typename detector_t::transform_container &transforms,
        typename detector_t::mask_container &masks,
        typename detector_t::geometry_context ctx = {}) const
        -> dindex_range = 0;
};

/// @brief Bind transform and mask data together for surface building.
///
/// Surface data is kept in separate detector containers and linked by indices,
/// so during geometry building, great care has to be taken to make sure that
/// all components of a surface get sorted and linked into the containers
/// together and associated with the correct surface.
template <typename detector_t, typename mask_shape_t, typename volume_link_t>
class surface_data_interface {
    public:
    surface_data_interface() = default;
    virtual ~surface_data_interface() = default;

    /// Uniform interface to surface data from different sources
    DETRAY_HOST
    virtual auto get_data()
        -> std::tuple<typename detector_t::transform3 &, volume_link_t &,
                      std::vector<typename detector_t::scalar_type> &> = 0;
};

}  // namespace detray