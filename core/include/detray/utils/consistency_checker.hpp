/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace detray::detail {

/// Checks every collection in a multistore to be emtpy and prints a warning
template <typename store_t, std::size_t... I>
void report_empty(const store_t &store, const std::string &store_name,
                  std::index_sequence<I...> /*seq*/) {

    ((store.template empty<store_t::value_types::to_id(I)>()
          ? std::cout << "WARNING: " << store_name
                      << " has empty collection no. " << I << std::endl
          : std::cout << ""),
     ...);
}

/// A functor that checks the surface descriptor and correct volume index in
/// every acceleration data structure for a given volume
struct surface_checker {

    /// Test the contained surfaces for consistency
    template <typename detector_t>
    DETRAY_HOST_DEVICE void operator()(
        const typename detector_t::surface_type &sf_descr,
        const detector_t &det, const dindex vol_idx) const {

        const auto sf = surface{det, sf_descr};
        std::stringstream err_stream{};

        if (not sf.self_check(err_stream)) {
            throw std::invalid_argument(err_stream.str());
        }

        if (sf.volume() != vol_idx) {
            err_stream << "ERROR: Incorrect volume index on surface: vol "
                       << vol_idx << ", sf: " << sf;

            throw std::invalid_argument(err_stream.str());
        }

        // Does the mask link to an existing volume?
        if (!detail::is_invalid_value(sf.volume_link()) &&
            (sf.volume_link() >= det.volumes().size())) {
            err_stream << "ERROR: Incorrect volume link to non-existent volume "
                       << sf.volume_link();
            throw std::invalid_argument(err_stream.str());
        }

        // Check that the same surface is registered in the detector surface
        // lookup
        const auto sf_from_lkp = surface{det, det.surface(sf.barcode())};
        if (not(sf_from_lkp == sf)) {
            err_stream << "ERROR: Surfaces in volume and detector lookups "
                       << "differ:\n In volume acceleration data structure: "
                       << sf << "\nIn detector surface lookup: " << sf_from_lkp;

            throw std::runtime_error(err_stream.str());
        }
    }

    /// Test wether a given surface @param check_desc is properly registered at
    /// least once in one of the volume acceleration data structures
    ///
    /// @param ref_descr one of the surfaces in the volumes acceleration data
    /// @param check_descr surface that we are searching for
    /// @param success communicate success to the outside
    template <typename detector_t>
    DETRAY_HOST_DEVICE void operator()(
        const typename detector_t::surface_type &ref_descr,
        const typename detector_t::surface_type &check_descr, bool &success,
        const detector_t &det) const {

        // Check that the surface is being searched for in the right volume
        // The volume index of the ref_descr must be checked to be correct
        // beforehand, e.g. by the call operator above
        if (ref_descr.volume() != check_descr.volume()) {
            std::stringstream err_stream{};
            err_stream << "Incorrect volume index on surface: "
                       << surface{det, check_descr};

            throw std::invalid_argument(err_stream.str());
        }

        // Check if it is the surface we are looking for
        if (ref_descr == check_descr) {
            success = true;
        }
    }
};

/// @brief Checks wether the data containers of a detector are empty
///
/// In case the default metadata is used, the unused containers are allowed to
/// be empty.
template <typename detector_t>
inline void check_empty(const detector_t &det) {

    // Check if there is at least one portal in the detector
    auto find_portals = [&det]() -> bool {
        if (det.portals().empty()) {
            return false;
        }
        // In the brute force finder, also other surfaces can be contained, e.g.
        // passive surfaces (depends on the detector)
        for (const auto &pt_desc : det.portals()) {
            if (pt_desc.is_portal()) {
                return true;
            }
        }
        return false;
    };

    // Check if there is at least one volume in the detector volume finder
    auto find_volumes =
        [](const typename detector_t::volume_finder &vf) -> bool {
        for (const auto &v : vf.all()) {
            if (not detail::is_invalid_value(v)) {
                return true;
            }
        }
        return false;
    };

    // Fatal errors
    if (det.volumes().empty()) {
        throw std::runtime_error("ERROR: No volumes in detector");
    }
    if (det.surfaces().empty()) {
        throw std::runtime_error("ERROR: No surfaces found");
    }
    if (det.transform_store().empty()) {
        throw std::runtime_error("ERROR: No transforms in detector");
    }
    if (det.mask_store().all_empty()) {
        throw std::runtime_error("ERROR: No masks in detector");
    }
    if (not find_portals()) {
        throw std::runtime_error("ERROR: No portals in detector");
    }

    // Warnings

    // Check for empty mask collections
    detail::report_empty(
        det.mask_store(), "mask store",
        std::make_index_sequence<detector_t::masks::n_types>{});

    // Check the material description
    if (det.material_store().all_empty()) {
        std::cout << "WARNING: No material in detector" << std::endl;
    } else {
        // Check for empty material collections
        detail::report_empty(
            det.material_store(), "material store",
            std::make_index_sequence<detector_t::materials::n_types>{});
    }

    // Check for empty acceleration data structure collections (e.g. grids)
    detail::report_empty(
        det.accelerator_store(), "acceleration data structures store",
        std::make_index_sequence<detector_t::accel::n_types>{});

    // Check volume search data structure
    if (not find_volumes(det.volume_search_grid())) {
        std::cout << "WARNING: No entries in volume finder" << std::endl;
    }
}

/// @brief Checks the internal consistency of a detector
template <typename detector_t>
inline bool check_consistency(const detector_t &det) {
    check_empty(det);

    std::stringstream err_stream{};
    // Check the volumes
    for (const auto &[idx, vol_desc] :
         detray::views::enumerate(det.volumes())) {
        const auto vol = detector_volume{det, vol_desc};

        // Check that nothing is obviously broken
        if (not vol.self_check(err_stream)) {
            throw std::invalid_argument(err_stream.str());
        }

        // Check consistency in the context of the owning detector
        if (vol.index() != idx) {
            err_stream << "ERROR: Incorrect volume index! Found volume:\n"
                       << vol << "\nat index " << idx;
            throw std::invalid_argument(err_stream.str());
        }

        // Go through the acceleration data structures and check the surfaces
        vol.template visit_surfaces<detail::surface_checker>(det, vol.index());
    }

    // Check the surfaces in the detector's surface lookup
    for (const auto &[idx, sf_desc] :
         detray::views::enumerate(det.surfaces())) {
        const auto sf = surface{det, sf_desc};

        // Check that nothing is obviously broken
        if (not sf.self_check(err_stream)) {
            err_stream << "\nat surface no. " << std::to_string(idx);
            throw std::invalid_argument(err_stream.str());
        }

        // Check consistency in the context of the owning detector
        if (sf.index() != idx) {
            err_stream << "ERROR: Incorrect surface index! Found surface:\n"
                       << sf << "\nat index " << idx;
            throw std::invalid_argument(err_stream.str());
        }

        // Check that the surface can be found in its volume's acceleration
        // data structures (if there are no grids, must at least be in the
        // brute force method)
        const auto vol = detector_volume{det, sf.volume()};
        bool is_registered = false;

        vol.template visit_surfaces<detail::surface_checker>(
            sf_desc, is_registered, det);

        if (not is_registered) {
            err_stream << "ERROR: Found surface that is not part of its "
                       << "volume's navigation acceleration data structures:\n"
                       << "Surface: " << sf;
            throw std::invalid_argument(err_stream.str());
        }
    }

    return true;
}

}  // namespace detray::detail
