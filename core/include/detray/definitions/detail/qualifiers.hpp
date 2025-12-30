/**
 * DETRAY library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if defined(__CUDACC__) || defined(__HIP__)
#define DETRAY_DEVICE __device__
#else
#define DETRAY_DEVICE
#endif

#if defined(__CUDACC__) || defined(__HIP__)
#define DETRAY_HOST __host__
#else
#define DETRAY_HOST
#endif

#if defined(__CUDACC__) || defined(__HIP__)
#define DETRAY_HOST_DEVICE __host__ __device__
#else
#define DETRAY_HOST_DEVICE
#endif

#if defined(__CUDACC__) || defined(__HIP__)
#define DETRAY_ALIGN(x) __align__(x)
#else
#define DETRAY_ALIGN(x) alignas(x)
#endif

#if defined(__CUDACC__) || defined(__HIP__) || defined(__GNUC__)
#define DETRAY_INLINE __attribute__((always_inline))
#else
#define DETRAY_INLINE
#endif

// @see
// https://stackoverflow.com/questions/78071873/gcc-preprocessor-macro-and-pragma-gcc-unroll
#if defined(__clang__)
#define ARG_TO_STRING(A) #A
#define DETRAY_UNROLL_N(n) _Pragma(ARG_TO_STRING(clang loop unroll_count(n)))
#elif defined(__GNUC__) || defined(__GNUG__)
#define ARG_TO_STRING(A) #A
#if __GNUC__ >= 14
#define DETRAY_UNROLL_N(n) _Pragma(ARG_TO_STRING(GCC unroll n))
#else
// For versions below 14, template parameters apparently cannot be used
#define DETRAY_UNROLL_N(n) _Pragma(ARG_TO_STRING(GCC unroll 8))
#endif
#else
// Unknown compiler or does not support unrolling directives
#define DETRAY_UNROLL_N(n)
#endif