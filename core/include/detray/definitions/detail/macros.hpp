// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

// Define a macro that assures there is no device compilation
#if not defined(__CUDACC__) && not defined(CL_SYCL_LANGUAGE_VERSION) && \
    not defined(SYCL_LANGUAGE_VERSION) && not defined(__HIP__)
#define __NO_DEVICE__
#endif
