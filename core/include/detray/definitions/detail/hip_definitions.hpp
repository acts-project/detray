/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if defined(__HIPCC__) || defined(__NVCC__)

#include <assert.h>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>  // for hipDeviceSynchronize, hipGetLastError
#include <stdio.h>
#include <stdlib.h>

/// Number of threads per Warp
#define WARP_SIZE 32

/// Helper macro used for checking  , type return values.

#define DETRAY_HIP_ERROR_CHECK(ans) \
    { hipAssert((ans), __FILE__, __LINE__); }
inline void hipAssert(hipError_t code, const char *file, int line,
                      bool abort = true) {
    if (code != hipSuccess) {
        fprintf(stderr, "HIPassert: %s %s %d\n", hipGetErrorString(code), file,
                line);
        if (abort)
            exit(code);
    }
}

#endif
