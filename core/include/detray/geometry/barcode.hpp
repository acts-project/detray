
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

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
    DETRAY_HOST_DEVICE
    constexpr explicit barcode(value_t encoded) : m_value(encoded) {}
    /// Construct default barcode with all values set to zero.
    constexpr barcode() = default;
    barcode(barcode&&) = default;
    barcode(const barcode&) = default;
    ~barcode() = default;
    barcode& operator=(barcode&&) = default;
    barcode& operator=(const barcode&) = default;

    /// Return the encoded value.
    DETRAY_HOST_DEVICE
    constexpr value_t value() const { return m_value; }

    /// Return the volume identifier.
    DETRAY_HOST_DEVICE
    constexpr value_t volume() const { return getBits(k_volume_mask); }

    /// Return the approach identifier.
    DETRAY_HOST_DEVICE
    constexpr surface_id id() const {
        return static_cast<surface_id>(getBits(k_id_mask));
    }

    /// Return the sensitive identifier.
    DETRAY_HOST_DEVICE
    constexpr value_t index() const { return getBits(k_index_mask); }

    /// Return the extra identifier
    /// Usage can be experiment-specific, like tagging which kind of detector a
    /// surface object corresponds to, or which subsystem it belongs to
    DETRAY_HOST_DEVICE
    constexpr value_t extra() const { return getBits(k_extra_mask); }

    /// Set the volume identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_volume(value_t volume) {
        return set_bits(k_volume_mask, volume);
    }

    /// Set the boundary identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_id(surface_id id) {
        return set_bits(k_id_mask, static_cast<value_t>(id));
    }

    /// Set the extra identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_index(value_t index) {
        return set_bits(k_index_mask, index);
    }

    /// Set the extra identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_extra(value_t extra) {
        return set_bits(k_extra_mask, extra);
    }

    private:
    // clang-format off
    static constexpr value_t k_volume_mask = 0xff00000000000000; // (2^8)-1 =
   255 volumes static constexpr value_t k_id_mask     = 0x00f0000000000000; //
   (2^4)-1 = 15 surface categories static constexpr value_t k_index_mask  =
   0x000fffffffffff00; // (2^44)-1 = 4095 surfaces static constexpr value_t
   k_extra_mask  = 0x00000000000000ff; // (2^8)-1 = 255 extra values
    // clang-format on

    value_t m_value{dindex_invalid};*/
// Use compiler builtin to extract the number of trailing bits from the
// mask. The builtin should be available on all supported compilers.
// Need unsigned long long version (...ll) to ensure uint64_t
// compatibility.
/*#if defined(__CUDACC__)
    /// Extract the bit shift necessary to access the masked values - CUDA
    ///
    /// @note undefined behaviour for mask == 0 which we should not have.
    DETRAY_DEVICE
    static constexpr int extractShift(value_t mask) {
        if (mask == k_volume_mask) { return 56; }
        else if (mask == k_id_mask) { return 52; }
        else if (mask == k_index_mask) { return 8; }
        else if (mask == k_extra_mask) { return 0; }
        else { assert(false && "invalid mask"); return -1; }
        // A '__ctz' builtin is not available in CUDA, so the mask needs to
        // be reversed and then the leading zeores are counted
        // return __clzll(__brevll(mask));
    }
#elif !defined(__CUDACC__)
    /// Extract the bit shift necessary to access the masked values -
host/SYCL
    ///
    /// @note undefined behaviour for mask == 0 which we should not have.
    DETRAY_HOST
    static constexpr int extractShift(value_t mask) {
        // should also be supported in SYCL:
https://github.com/intel/llvm/blob/sycl/clang/docs/LanguageExtensions.rst#id91
        return __builtin_ctzll(mask);
    }
#endif*/

/// Extract the bit shift necessary to access the masked values.
/*DETRAY_HOST_DEVICE
static constexpr int extractShift(value_t mask) {
    // use compiler builtin to extract the number of trailing bits from the
    // mask. the builtin should be available on all supported compilers.
    // need unsigned long long version (...ll) to ensure uint64_t
    // compatibility.
    // WARNING undefined behaviour for mask == 0 which we should not have.
    return __builtin_ctzll(mask);
}

/// Extract the masked bits from the encoded value.
DETRAY_HOST_DEVICE
constexpr value_t getBits(value_t mask) const {
    return (m_value & mask) >> extractShift(mask);
}

/// Set the masked bits to id in the encoded value.
DETRAY_HOST_DEVICE
constexpr barcode& set_bits(value_t mask, value_t id) {
    m_value = (m_value & ~mask) | ((id << extractShift(mask)) & mask);
    // return *this here so we need to write less lines in the 'set' methods
    return *this;
}

DETRAY_HOST_DEVICE
friend constexpr bool operator==(barcode lhs, barcode rhs) {
    return lhs.m_value == rhs.m_value;
}

DETRAY_HOST_DEVICE
friend constexpr bool operator!=(barcode lhs, barcode rhs) {
    return lhs.m_value != rhs.m_value;
}

DETRAY_HOST_DEVICE
friend constexpr bool operator<(barcode lhs, barcode rhs) {
    return lhs.m_value < rhs.m_value;
}

DETRAY_HOST
friend std::ostream& operator<<(std::ostream& os, const barcode c) {
    // zero represents an invalid/undefined identifier
    if (c.volume() == dindex_invalid) {
        return (os << "undefined");
    }

    static const char* const names[] = {
        "vol = ", "id = ", "index = ", "extra = "};
    const geometry::barcode::value_t levels[] = {
        c.volume(), static_cast<dindex>(c.id()), c.index(), c.extra()};

    bool writeSeparator = false;
    for (auto i = 0u; i < 4u; ++i) {
        if (writeSeparator) {
            os << " | ";
        }
        os << names[i] << levels[i];
        writeSeparator = true;
    }
    return os;
}
};

}  // namespace detray::geometry

namespace std {

// specialize std::hash so barcode can be used e.g. in an unordered_map
template <>
struct hash<detray::geometry::barcode> {
auto operator()(detray::geometry::barcode gid) const noexcept {
    return std::hash<detray::geometry::barcode::value_t>()(gid.value());
}
};

}  // namespace std*/

class barcode {

    public:
    using value_t = dindex;

    /// Construct from an already encoded value.
    // constexpr explicit barcode(value_t encoded) : m_value(encoded) {}
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
    /// Return the sensitive identifier.
    DETRAY_HOST_DEVICE
    constexpr value_t extra() const { return m_extra; }

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
    /// Set the extra identifier.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_extra(value_t extra) {
        m_extra = extra;
        return *this;
    }

    /// Check wether the barcode is valid to use.
    DETRAY_HOST_DEVICE
    constexpr bool is_invalid() const {
        return ((m_volume == 255UL) | (m_id == 15UL) |
                (m_index == std::pow(2UL, 44UL) - 1UL));
    }

    private:
    // clang-format off
    static constexpr value_t k_volume_mask = 0xff00000000000000; // (2^8)-1 = 255 volumes 
    static constexpr value_t k_id_mask     = 0x00f0000000000000; // (2^4)-1 = 15 surface categories 
    static constexpr value_t k_index_mask  = 0x000fffffffffff00; // (2^44)-1 = 4095 surfaces 
    static constexpr value_t k_extra_mask  = 0x00000000000000ff; // (2^8)-1 = 255 extra values
    // clang-format on

    value_t m_volume = 255UL, m_id = 15UL, m_index = std::pow(2UL, 44UL) - 1UL,
            m_extra = 255UL;

    value_t m_value{dindex_invalid};

    /// Extract the bit shift necessary to access the masked values.
    DETRAY_HOST_DEVICE
    static constexpr int extractShift(value_t mask) {
        // use compiler builtin to extract the number of trailing bits from the
        // mask. the builtin should be available on all supported compilers.
        // need unsigned long long version (...ll) to ensure uint64_t
        // compatibility.
        // WARNING undefined behaviour for mask == 0 which we should not have.
        return __builtin_ctzll(mask);
    }

    /// Extract the masked bits from the encoded value.
    DETRAY_HOST_DEVICE
    constexpr value_t getBits(value_t mask) const {
        return (m_value & mask) >> extractShift(mask);
    }

    /// Set the masked bits to id in the encoded value.
    DETRAY_HOST_DEVICE
    constexpr barcode& set_bits(value_t mask, value_t id) {
        // return *this here so we need to write less lines in the 'set' methods
        return *this;
    }

    DETRAY_HOST_DEVICE
    friend constexpr bool operator==(barcode lhs, barcode rhs) {
        return lhs.m_volume == rhs.m_volume and lhs.m_id == rhs.m_id and
               lhs.m_index == rhs.m_index and lhs.m_extra == rhs.m_extra;
    }
    DETRAY_HOST_DEVICE
    friend constexpr bool operator!=(barcode lhs, barcode rhs) {
        return lhs.m_volume != rhs.m_volume and lhs.m_id != rhs.m_id and
               lhs.m_index != rhs.m_index and lhs.m_extra != rhs.m_extra;
    }
    DETRAY_HOST_DEVICE
    friend constexpr bool operator<(barcode lhs, barcode rhs) {
        return lhs.m_index < rhs.m_index;
    }
    DETRAY_HOST
    friend std::ostream& operator<<(std::ostream& os, const barcode c) {
        // zero represents an invalid/undefined identifier
        if (c.is_invalid()) {
            return (os << "undefined");
        }

        static const char* const names[] = {
            "vol = ", "id = ", "index = ", "extra = "};
        const geometry::barcode::value_t levels[] = {
            c.volume(), static_cast<dindex>(c.id()), c.index(), c.extra()};

        bool writeSeparator = false;
        for (auto i = 0u; i < 4u; ++i) {
            if (writeSeparator) {
                os << " | ";
            }
            os << names[i] << levels[i];
            writeSeparator = true;
        }
        return os;
    }
};

}  // namespace detray::geometry