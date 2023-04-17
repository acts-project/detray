/** Detray library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if __cpp_concepts >= 201907L
#define CONSTRAINT(x) x
#define DETRAY_HAVE_CONCEPTS
#else
#define CONSTRAINT(x) typename
#endif
