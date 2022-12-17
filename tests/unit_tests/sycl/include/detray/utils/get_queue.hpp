/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "queue_wrapper.hpp"

// SYCL include(s).
#include <CL/sycl.hpp>

namespace detray::sycl::details {

/// Helper function for getting a @c sycl::queue out of
/// @c detray::sycl::queue_wrapper (non-const)
::sycl::queue& get_queue(detray::sycl::queue_wrapper& queue) {
    assert(queue.queue() != nullptr);
    return *(reinterpret_cast<::sycl::queue*>(queue.queue()));
}

/// Helper function for getting a @c sycl::queue out of
/// @c detray::sycl::queue_wrapper (const)
const ::sycl::queue& get_queue(const detray::sycl::queue_wrapper& queue) {
    assert(queue.queue() != nullptr);
    return *(reinterpret_cast<const ::sycl::queue*>(queue.queue()));
}

}  // namespace detray::sycl::details
