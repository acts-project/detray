
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"

// System include(s)
#include <cstdint>
#include <iosfwd>
#include <ostream>
#include <utility>

namespace detray::geometry {

/// @brief Unique identifier for geometry objects, based on the ACTS
/// GeometryIdentifier
///
/// Encodes the volume index, the type of surface (portal, sensitive, passive
/// etc), an index to identify a surface in a geometry accelerator structure,
/// as well as two extra bytes that can be used to tag surfaces.
///
/// @note the detray barcode is not compatible with the ACTS GeometryIdentifier
///
/// @see
/// https://github.com/acts-project/acts/blob/main/Core/include/Acts/Geometry/GeometryIdentifier.hpp
/*class barcode {

    public:
    using value_t = uint64_t;

    /// Construct from an already encoded value.
    constexpr explicit barcode(value_t encoded) : m_value(encoded) {}
    /// Construct default barcode with all values set to zero.
    constexpr barcode() = default;
    barcode(barcode&&) = default;
    barcode(const barcode&) = default;
    ~barcode() = default;
    barcode& operator=(barcode&&) = default;
    barcode& operator=(const barcode&) = default;

    /// Return the encoded value.
    constexpr value_t value() const { return m_value; }

    /// Return the volume identifier.
    constexpr value_t volume() const { return getBits(k_volume_mask); }
    /// Return the approach identifier.
    constexpr value_t id() const { return getBits(k_id_mask); }
    /// Return the sensitive identifier.
    constexpr value_t index() const { return getBits(k_index_mask); }
    /// Return the extra identifier
    /// Usage can be experiment-specific, like tagging which kind of detector a
    /// surface object corresponds to, or which subsystem it belongs to
    constexpr value_t extra() const { return getBits(k_extra_mask); }

    /// Set the volume identifier.
    constexpr barcode& set_volume(value_t volume) {
        return set_bits(k_volume_mask, volume);
    }
    /// Set the boundary identifier.
    constexpr barcode& set_id(value_t id) {
        return set_bits(k_id_mask, id);
    }
    /// Set the extra identifier
    constexpr barcode& set_index(value_t index) {
        return set_bits(k_index_mask, index);
    }
    /// Set the extra identifier
    constexpr barcode& set_extra(value_t extra) {
        return set_bits(k_extra_mask, extra);
    }

    private:
    // clang-format off
    static constexpr value_t k_volume_mask = 0xff00000000000000; // (2^8)-1 =
255 volumes static constexpr value_t k_id_mask     = 0x00f0000000000000; //
(2^4)-1 = 16 surface categories static constexpr value_t k_index_mask  =
0x000fffffffffff00; // (2^44)-1 = 4095 surfaces static constexpr value_t
k_extra_mask  = 0x00000000000000ff; // (2^8)-1 = 255 extra values
    // clang-format on

    value_t m_value = 0;

    /// Extract the bit shift necessary to access the masked values.
    static constexpr int extractShift(value_t mask) {
        // use compiler builtin to extract the number of trailing bits from the
        // mask. the builtin should be available on all supported compilers.
        // need unsigned long long version (...ll) to ensure uint64_t
        // compatibility.
        // WARNING undefined behaviour for mask == 0 which we should not have.
        return __builtin_ctzll(mask);
    }
    /// Extract the masked bits from the encoded value.
    constexpr value_t getBits(value_t mask) const {
        return (m_value & mask) >> extractShift(mask);
    }
    /// Set the masked bits to id in the encoded value.
    constexpr barcode& set_bits(value_t mask, value_t id) {
        m_value = (m_value & ~mask) | ((id << extractShift(mask)) & mask);
        // return *this here so we need to write less lines in the set...
methods return *this;
    }

    friend constexpr bool operator==(barcode lhs,
                                    barcode rhs) {
        return lhs.m_value == rhs.m_value;
    }
    friend constexpr bool operator!=(barcode lhs,
                                    barcode rhs) {
        return lhs.m_value != rhs.m_value;
    }
    friend constexpr bool operator<(barcode lhs,
                                    barcode rhs) {
        return lhs.m_value < rhs.m_value;
    }
};

std::ostream& operator<<(std::ostream& os, barcode c) {
  // zero represents an invalid/undefined identifier
  if (c.value() == dindex_invalid) {
    return (os << "undefined");
  }

  static const char* const names[] = {
      "vol=", "id=", "index=", "ext=",
  };
  const barcode::value_t levels[] = {c.volume(), c.id(), c.index(), id.extra()};

  bool writeSeparator = false;
  for (const auto i = 0u; i < 5u; ++i) {
    if (levels[i] != 0u) {
      if (writeSeparator) {
        os << '|';
      }
      os << names[i] << levels[i];
      writeSeparator = true;
    }
  }
  return os;
}

}  // namespace detray::geometry


namespace std {

// specialize std::hash so barcode can be used e.g. in an unordered_map
template <>
struct hash<detray::geometry::barcode> {
    auto operator()(detray::geometry::barcode gid) const noexcept {
        return std::hash<detray::geometry::barcode::value_t>()(gid.value());
    }
};

} // namespace std*/

class barcode {

    public:
    using value_t = dindex;

    /// Construct from an already encoded value.
    constexpr explicit barcode(value_t encoded) : m_index(encoded) {}
    /// Construct default barcode with all values set to zero.
    constexpr barcode() = default;
    barcode(barcode&&) = default;
    barcode(const barcode&) = default;
    ~barcode() = default;
    barcode& operator=(barcode&&) = default;
    barcode& operator=(const barcode&) = default;

    /// Return the volume identifier.
    DETRAY_HOST_DEVICE
    constexpr value_t volume() const { return m_volume; }
    /// Return the approach identifier.
    DETRAY_HOST_DEVICE
    constexpr surface_id id() const { return static_cast<surface_id>(m_id); }
    /// Return the sensitive identifier.
    DETRAY_HOST_DEVICE
    constexpr value_t index() const { return m_index; }

    /// Set the volume identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_volume(value_t volume) {
        m_volume = volume;
        return *this;
    }
    /// Set the boundary identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_id(surface_id id) {
        m_id = static_cast<dindex>(id);
        return *this;
    }
    /// Set the extra identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_index(value_t index) {
        m_index = index;
        return *this;
    }

    private:
    value_t m_volume = dindex_invalid, m_id = dindex_invalid,
            m_index = dindex_invalid;

    DETRAY_HOST_DEVICE
    friend constexpr bool operator==(barcode lhs, barcode rhs) {
        return lhs.m_volume == rhs.m_volume and lhs.m_id == rhs.m_id and
               lhs.m_index == rhs.m_index;
    }
    DETRAY_HOST_DEVICE
    friend constexpr bool operator!=(barcode lhs, barcode rhs) {
        return lhs.m_volume != rhs.m_volume and lhs.m_id != rhs.m_id and
               lhs.m_index != rhs.m_index;
    }
    DETRAY_HOST_DEVICE
    friend constexpr bool operator<(barcode lhs, barcode rhs) {
        return lhs.m_index < rhs.m_index;
    }
    DETRAY_HOST
    friend std::ostream& operator<<(std::ostream& os, const barcode c) {
        // zero represents an invalid/undefined identifier
        if (c.volume() == dindex_invalid) {
            return (os << "undefined");
        }

        static const char* const names[] = {"vol=", "id=", "index="};
        const geometry::barcode::value_t levels[] = {
            c.volume(), static_cast<dindex>(c.id()), c.index()};

        bool writeSeparator = false;
        for (auto i = 0u; i < 3u; ++i) {
            if (levels[i] != 0u) {
                if (writeSeparator) {
                    os << '|';
                }
                os << names[i] << levels[i];
                writeSeparator = true;
            }
        }
        return os;
    }
};

}  // namespace detray::geometry