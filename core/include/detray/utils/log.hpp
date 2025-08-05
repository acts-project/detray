/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#ifdef DETRAY_ENABLE_LOGGING
#include <cxxabi.h>

#include <iostream>
#include <regex>
#include <string>

namespace detray::log::detail {
template <typename T>
inline std::string_view cpp_demangle() {
    static const std::string type_name = [] {
        const char *abiName = typeid(T).name();
        int status = 0;
        char *ret = abi::__cxa_demangle(abiName, nullptr, nullptr, &status);

        std::string s{ret};
        free(static_cast<void *>(ret));

        std::regex re{"detray::"};
        s = std::regex_replace(s, re, "");

        // Special case type-list for readability
        // Only attempt if the full type is a list
        std::regex re_list{R"(^types::list<(.*)>$)"};

        std::smatch match;
        if (std::regex_match(s, match, re_list)) {
            s = "[" + std::string{match[1]} + "]";
        } else {
        }
        return s;
    }();

    return type_name;
}
}  // namespace detray::log::detail

#define DETRAY_TYPENAME(type) detray::log::detail::cpp_demangle<type>()

#define __FILENAME__ \
    (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define DETRAY_LOG(lvl, x)                                                  \
    std::cout << __FILENAME__ << ":" << __LINE__ << " " << lvl << ": " << x \
              << std::endl;

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

#else
#define DETRAY_LOG(x)
#define DETRAY_TYPENAME(type)
#endif

#define DETRAY_DEBUG(x) DETRAY_LOG("DEBUG", x)
#define DETRAY_ERROR(x) DETRAY_LOG("ERROR", x)
