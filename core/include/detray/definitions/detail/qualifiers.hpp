// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

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
