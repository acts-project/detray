/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/** Helper utility to expand gens<N> into seq<0,1,2,...,N-1>
 *
 * Reference:
 * https://stackoverflow.com/questions/36612596/tuple-to-parameter-pack
 *
 * Usage with parameter pack of (Args...) with N parameters:
 *
 * typename gens<sizeof...(Args)>::type() => seq<0,1,2,...,N-1>
 */

template <unsigned int...>
struct seq {};

template <unsigned int N, unsigned int... S>
struct gens : gens<N - 1, N - 1, S...> {};

template <unsigned int... S>
struct gens<0, S...> {
    typedef seq<S...> type;
};

}  // namespace detray
