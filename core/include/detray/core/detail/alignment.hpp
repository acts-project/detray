// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

namespace detray::detail {

/// Creates detector view using "static" detector components and
/// a "misaligned" transform store
template <typename host_detector_type>
typename host_detector_type::view_type misaligned_detector_view(
    typename host_detector_type::buffer_type& det_buffer,
    typename host_detector_type::transform_container::buffer_type& trf_buffer) {
    typename host_detector_type::view_type detview{
        detray::get_data(
            detray::detail::get<0>(det_buffer.m_buffer)),  // volumes
        detray::get_data(
            detray::detail::get<1>(det_buffer.m_buffer)),  // surfaces
        detray::get_data(trf_buffer),                      // transforms
        detray::get_data(detray::detail::get<3>(det_buffer.m_buffer)),  // masks
        detray::get_data(
            detray::detail::get<4>(det_buffer.m_buffer)),  // materials
        detray::get_data(
            detray::detail::get<5>(det_buffer.m_buffer)),  // accelerators
        detray::get_data(detray::detail::get<6>(
            det_buffer.m_buffer))};  // volume search grid
    return detview;
}

}  // namespace detray::detail
