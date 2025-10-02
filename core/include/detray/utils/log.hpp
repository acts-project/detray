/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if DETRAY_LOG_LVL >= 0
// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <cstring>
#include <iostream>
#include <regex>
#include <source_location>
#include <string>

namespace detray::log::detail {

/// @returns the name of the current source file without the full path
/// @see
/// https://stackoverflow.com/questions/31050113/how-to-extract-the-source-filename-without-path-and-suffix-at-compile-time
constexpr const char *source_file_name(const char *path) {
    const char *file = path;
    while (*path) {
        if (*path++ == '/') {
            file = path;
        }
    }
    return file;
}

/// Print a type name for type @tparam T in the logs
template <typename T>
inline std::string_view process_typename() {
    static const std::string type_name = [] {
        std::string s{""};
        try {
            s = detray::types::demangle_type_name<T>();
        } catch (...) {
            return std::string{"unknown"};
        }

        if (s.empty()) {
            return std::string{"unknown"};
        }

        std::regex re{"detray::"};
        s = std::regex_replace(s, re, "");

        // Special case type-list for readability
        // Only attempt if the full type is a list
        std::regex re_list{R"(^types::list<(.*)>$)"};

        if (std::smatch match; std::regex_match(s, match, re_list)) {
            s = "[" + std::string{match[1]} + "]";
        }

        return s;
    }();

    return type_name;
}

}  // namespace detray::log::detail

// Define a macro that enables device logging using printf
#if defined(__CUDACC__) || defined(__HIP__)
#define __DEVICE_LOGGING__
#endif

#ifdef __DEVICE_LOGGING__
#define DETRAY_TYPENAME(type) "unknown type"
#else
#define DETRAY_TYPENAME(type) detray::log::detail::process_typename<type>()
#endif

#ifdef __DEVICE_LOGGING__
#define DETRAY_LOG_VECTOR(x) ""
#else
#define DETRAY_LOG_VECTOR(x)              \
    [&]() {                               \
        std::stringstream _vec_os;        \
        std::size_t _vec_i = 0;           \
        for (const auto &_vec_elem : x) { \
            if (_vec_i > 0) {             \
                _vec_os << ", ";          \
            }                             \
            _vec_os << _vec_elem;         \
            _vec_i++;                     \
        }                                 \
        return _vec_os.str();             \
    }()
#endif

// String that represents the backend that emitted the log message
#if defined(__CUDACC__)
#define __BACKEND__ "CUDA"
#elif defined(__HIP__)
#define __BACKEND__ "HIP"
#elif defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#define __BACKEND__ "SYCL"
#else
#define __BACKEND__ "HOST"
#endif

#ifdef __DEVICE_LOGGING__
#define __FILENAME__ detray::log::detail::source_file_name(__FILE__)
#else
#define __FILENAME__ \
    (std::strrchr(__FILE__, '/') ? std::strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

// Print 'x' in the host logs only
#ifndef __DEVICE_LOGGING__
#define DETRAY_LOG_STREAM(lvl, x)                                       \
    std::clog << "DETRAY " << lvl << " (HOST): " << __FILENAME__ << ":" \
              << std::source_location::current().line() << " " << x     \
              << std::endl
#else
#define DETRAY_LOG_STREAM(lvl, x)
#endif

// Print 'x' in printf
// @note 'x' is the format string and the variadic arguments the corresponding
// values
// @note logging is currently disabled for SYCL builds
#if !defined(CL_SYCL_LANGUAGE_VERSION) && !defined(SYCL_LANGUAGE_VERSION)
#ifdef __DEVICE_LOGGING__
#define DETRAY_LOG_PRINTF(lib, lvl, x, ...)                                  \
    printf("%s %s (%s): %s:%d " x "\n", lib, lvl, __BACKEND__, __FILENAME__, \
           __LINE__ __VA_OPT__(, ) __VA_ARGS__);
#else
#define DETRAY_LOG_PRINTF(lib, lvl, x, ...)                             \
    fprintf(stderr, "%s %s (%s): %s:%d " x "\n", lib, lvl, __BACKEND__, \
            __FILENAME__, __LINE__ __VA_OPT__(, ) __VA_ARGS__);
#endif
#else  // ifndef SYCL
#define DETRAY_LOG_PRINTF(lib, lvl, x, ...)
#endif

#else  // DETRAY_LOG_LVL < 0
#define DETRAY_LOG_STREAM(lvl, x)
#define DETRAY_LOG_PRINTF(lib, lvl, x, ...)
#define DETRAY_TYPENAME(type)
#endif

// Warnings
#define DETRAY_WARN_HOST(x) DETRAY_LOG_STREAM("WARNING", x)
#ifdef __DEVICE_LOGGING__
#define DETRAY_WARN_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "WARNING", x, __VA_ARGS__)
#else
#define DETRAY_WARN_DEVICE(x, ...)
#endif
#define DETRAY_WARN_HOST_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "WARNING", x, __VA_ARGS__)

#define DETRAY_ERROR_HOST(x) DETRAY_LOG_STREAM("ERROR", x)
#ifdef __DEVICE_LOGGING__
#define DETRAY_ERROR_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "ERROR", x, __VA_ARGS__)
#else
#define DETRAY_ERROR_DEVICE(x, ...)
#endif
#define DETRAY_ERROR_HOST_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "ERROR", x, __VA_ARGS__)

#define DETRAY_FATAL_HOST(x) DETRAY_LOG_STREAM("FATAL", x)
#ifdef __DEVICE_LOGGING__
#define DETRAY_FATAL_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "FATAL", x, __VA_ARGS__)
#else
#define DETRAY_FATAL_DEVICE(x, ...)
#endif
#define DETRAY_FATAL_HOST_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "FATAL", x, __VA_ARGS__)

// Info
#if DETRAY_LOG_LVL > 0
#define DETRAY_INFO_HOST(x) DETRAY_LOG_STREAM("INFO", x)
#ifdef __DEVICE_LOGGING__
#define DETRAY_INFO_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "INFO", x, __VA_ARGS__)
#else
#define DETRAY_INFO_DEVICE(x, ...)
#endif
#define DETRAY_INFO_HOST_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "INFO", x, __VA_ARGS__)
#else
#define DETRAY_INFO_HOST(x)
#define DETRAY_INFO_DEVICE(x, ...)
#define DETRAY_INFO_HOST_DEVICE(x, ...)
#endif

// Info
#if DETRAY_LOG_LVL > 1
#define DETRAY_VERBOSE_HOST(x) DETRAY_LOG_STREAM("VERBOSE", x)
#ifdef __DEVICE_LOGGING__
#define DETRAY_VERBOSE_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "VERBOSE", x, __VA_ARGS__)
#else
#define DETRAY_VERBOSE_DEVICE(x, ...)
#endif
#define DETRAY_VERBOSE_HOST_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "VERBOSE", x, __VA_ARGS__)
#else
#define DETRAY_VERBOSE_HOST(x)
#define DETRAY_VERBOSE_DEVICE(x, ...)
#define DETRAY_VERBOSE_HOST_DEVICE(x, ...)
#endif

// Debug
#if DETRAY_LOG_LVL > 2
#define DETRAY_DEBUG_HOST(x) DETRAY_LOG_STREAM("DEBUG", x)
#ifdef __DEVICE_LOGGING__
#define DETRAY_DEBUG_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "DEBUG", x, __VA_ARGS__)
#else
#define DETRAY_DEBUG_DEVICE(x, ...)
#endif
#define DETRAY_DEBUG_HOST_DEVICE(x, ...) \
    DETRAY_LOG_PRINTF("DETRAY", "DEBUG", x, __VA_ARGS__)
#else
#define DETRAY_DEBUG_HOST(x)
#define DETRAY_DEBUG_DEVICE(x, ...)
#define DETRAY_DEBUG_HOST_DEVICE(x, ...)
#endif
