/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

// System include(s).
#include <cstddef>

namespace detray::algebra::fastor {

/// @brief Fastor vector type that provides equality op
template <concepts::scalar T, std::size_t N>
class Vector : public Fastor::Tensor<T, N> {

    public:
    /// Inherit all constructors from the base class
    using Fastor::Tensor<T, N>::Tensor;

    private:
    /// Equality operator for fastor vectors
    template <std::size_t M, concepts::scalar S>
    constexpr friend bool operator==(const Vector<S, M>& lhs,
                                     const Vector<S, M>& rhs);
    /// @}

};  // class Vector

template <std::size_t M, concepts::scalar S>
constexpr bool operator==(const Vector<S, M>& lhs, const Vector<S, M>& rhs) {
    return Fastor::isequal(static_cast<Fastor::Tensor<S, M>>(lhs),
                           static_cast<Fastor::Tensor<S, M>>(rhs), 0.00001f);
}

}  // namespace detray::algebra::fastor
