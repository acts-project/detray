/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <any>
#include <utility>

namespace detray {

namespace detail {

/// @brief ACTS-style context object.
///
/// @see
/// https://github.com/acts-project/acts/blob/main/Core/include/Acts/Utilities/detail/ContextType.hpp
class data_context {
    /*public:

   /// Does nothing.
   constexpr data_context() noexcept {};

   /// Move construct a new Context Type object from anything. Must be explicit.
   ///
   /// @tparam T The type of the value to construct from
   /// @param value The value to construct from
   template <typename T>
   DETRAY_HOST_DEVICE
   explicit data_context(T&& value) : m_data{std::move(value)} {}

   /// Copy construct a new Context Type object from anything. Must be explicit.
   ///
   /// @tparam T The type of the value to construct from
   /// @param value The value to construct from
   template <typename T>
   DETRAY_HOST_DEVICE
   explicit data_context(const T& value) : m_data{value} {}

   template<class T, class... Args>
   DETRAY_HOST_DEVICE
   explicit data_context(std::in_place_type_t<T>, Args&&... args)
       : m_data{std::in_place, std::forward<Args>(args)...} {}

   /// Move assignment of anything to this object is allowed.
   ///
   /// @tparam T The type of the value to assign
   /// @param value The value to assign
   /// @return data_context&
   template <typename T>
   DETRAY_HOST_DEVICE
   auto operator=(T&& value) noexcept -> data_context&  {
       m_data = std::move(value);
       return *this;
   }

   /// Copy assignment of anything to this object is allowed.
   ///
   /// @tparam T The type of the value to assign
   /// @param value The value to assign
   /// @return data_context&
   template <typename T>
   DETRAY_HOST_DEVICE
   auto operator=(const T& value) -> data_context& {
       m_data = value;
       return *this;
   }

   /// Retrieve a reference to the contained type
   ///
   /// @tparam T The type to attempt to retrieve the value as
   /// @return Reference to the contained value
   template <typename T>
   DETRAY_HOST_DEVICE
   auto get() -> std::decay_t<T>& {
       return std::any_cast<std::decay_t<T>&>(m_data);
   }

   /// Retrieve a reference to the contained type
   ///
   /// @tparam T The type to attempt to retrieve the value as
   /// @return Reference to the contained value
   template <typename T>
   DETRAY_HOST_DEVICE
   auto get() const -> const std::decay_t<T>&  {
       return std::any_cast<const std::decay_t<T>&>(m_data);
   }

   /// Check if the contained type is initialized.
   /// @return Boolean indicating whether a type is present
   DETRAY_HOST_DEVICE
   auto has_value() const noexcept -> bool { return m_data.has_value(); }

   private:
   std::any m_data;*/
};

}  // namespace detail

/// Placeholder context type
class empty_context {};

/// Context type for geometry data
class geometry_context : public detail::data_context {};

/// Context type for magnetic field data
struct magnetic_field_context : public detail::data_context {};

}  // namespace detray
