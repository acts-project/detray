/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/surface_finders/grid/axis.hpp"

namespace detray {

/// @brief Transforms a local axis bin index to a global grid bin index and vice
/// versa
///
/// Serializers allow to create a memory local data layout if advantegeous.
///
/// @note the bin indices are expected to start at zero.
template <std::size_t kDIM>
struct simple_serializer {};

/// @brief Simple serializer specialization for a single axis
template <>
struct simple_serializer<1> {

    /// @returns the axis local bin, which is also the global bin
    template <typename multi_axis_t>
    DETRAY_HOST_DEVICE auto operator()(multi_axis_t & /*axes*/,
                                       n_axis::multi_bin<1> mbin) const
        -> dindex {
        return mbin[0];
    }

    /// @returns the global bin, which is also the axis local bin
    template <typename multi_axis_t,
              template <typename, std::size_t> class array_t = darray>
    DETRAY_HOST_DEVICE auto operator()(multi_axis_t & /*axes*/,
                                       dindex gbin) const
        -> n_axis::multi_bin<1> {
        return {gbin};
    }
};

/// @brief Simple serializer specialization for a 3D multi-axis
template <>
struct simple_serializer<2> {

    /// @brief Create a serial bin from a multi-bin - 3D
    ///
    /// @tparam multi_axis_t is the type of multi-dimensional axis
    ///
    /// @param axes contains all axes (multi-axis)
    /// @param mbin contains a bin index for every axis in the multi-axis.
    ///
    /// @returns a dindex for the bin data storage
    template <typename multi_axis_t>
    DETRAY_HOST_DEVICE auto operator()(multi_axis_t &axes,
                                       n_axis::multi_bin<2> mbin) const
        -> dindex {
        dindex offset{mbin[1] * axes.template get_axis<0>().nbins()};
        return offset + mbin[0];
    }

    /// @brief Create a bin tuple from a serialized bin - 3D
    ///
    /// @tparam multi_axis_t is the type of multi-dimensional axis
    ///
    /// @param axes contains all axes (multi-axis)
    /// @param gbin the global (serial) bin
    ///
    /// @return a 2-dimensional multi-bin
    template <typename multi_axis_t>
    DETRAY_HOST_DEVICE auto operator()(multi_axis_t &axes, dindex gbin) const
        -> n_axis::multi_bin<2> {
        dindex nbins_axis0 = axes.template get_axis<0>().nbins();

        dindex bin0{gbin % nbins_axis0};
        dindex bin1{static_cast<dindex>(gbin / nbins_axis0)};

        return {bin0, bin1};
    }
};

/// @brief Simple serializer specialization for a 3D multi-axis
template <>
struct simple_serializer<3> {

    /// @brief Create a serial bin from a multi-bin - 3D
    ///
    /// @tparam multi_axis_t is the type of multi-dimensional axis
    ///
    /// @param axes contains all axes (multi-axis)
    /// @param mbin contains a bin index for every axis in the multi-axis.
    ///
    /// @returns a dindex for the bin data storage
    template <typename multi_axis_t>
    DETRAY_HOST_DEVICE auto operator()(multi_axis_t &axes,
                                       n_axis::multi_bin<3> mbin) const
        -> dindex {
        dindex nbins_axis0 = axes.template get_axis<0>().nbins();
        dindex nbins_axis1 = axes.template get_axis<1>().nbins();

        dindex offset1{mbin[1] * nbins_axis0};
        dindex offset2{mbin[2] * (nbins_axis0 * nbins_axis1)};

        return offset1 + offset2 + mbin[0];
    }

    /// @brief Create a bin tuple from a serialized bin - 3D
    ///
    /// @tparam multi_axis_t is the type of multi-dimensional axis
    ///
    /// @param axes contains all axes (multi-axis)
    /// @param gbin the global (serial) bin
    ///
    /// @return a 3-dimensional multi-bin
    template <typename multi_axis_t>
    DETRAY_HOST_DEVICE auto operator()(multi_axis_t &axes, dindex gbin) const
        -> n_axis::multi_bin<3> {
        dindex nbins_axis0 = axes.template get_axis<0>().nbins();
        dindex nbins_axis1 = axes.template get_axis<1>().nbins();

        dindex bin0{gbin % nbins_axis0};
        dindex bin1{static_cast<dindex>(gbin / (nbins_axis0)) % nbins_axis1};
        dindex bin2{static_cast<dindex>(gbin / (nbins_axis0 * nbins_axis1))};
        return {bin0, bin1, bin2};
    }
};

}  // namespace detray
