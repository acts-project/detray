# detray

Detray is part of the ACTS project (R&D line for parallelization), the ACTS project can be found: https://github.com/acts-project/acts.

This is a C++17 header only library for detector surface intersections using different algebra plugin libraries. It follows the navigation and propagation concept of ACTS, however, with an attempt to create
a geometry without polymorphic inheritance structure.

## Repository structure

The code for defining detector surfaces is within the `core` directory, the algebra plugins are in the `plugin` directory, and tests are all within the `test` library.

Testing and benchmarking is done with `googletest` and `google/benchmark`.


## Getting started

The respository is meant to be possible to build "out of the box", with standard
CMake build procedures.

```shell
git clone https://github.com/acts-project/detray.git
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -S detray -B detray-build
cmake --build detray-build
```

The following cache variables are available to influence which parts of the
project would be built:

- `DETRAY_EIGEN_PLUGIN` (`ON` by default), `DETRAY_SMATRIX_PLUGIN`
  (`OFF` by default), `DETRAY_VC_PLUGIN` (`ON` by default): Boolean
  flags turning the build of [Eigen](https://eigen.tuxfamily.org),
  [SMatrix](https://root.cern/doc/master/group__SMatrixGroup.html) and
  [Vc](https://github.com/VcDevel/Vc) using code on or off.
  * Note that [Algebra Plugins](https://github.com/acts-project/algebra-plugins)
    must have all of the appropriate options enabled for whichever option
    is turned on from these.
- `DETRAY_DISPLAY`: Boolean option turning on the build of `detray::display`,
  and additional helpers for displaying a geometry (`OFF` by default);
- `DETRAY_BUILD_CUDA`: Boolean option turning on the build of all CUDA code
  in the repository (`ON` by default, if CUDA is available);
- `DETRAY_BUILD_TESTING`: Turn the build/setup of the unit tests on/off
  (`ON` by default);
  * `DETRAY_BENCHMARKS`: Boolean option turning on the build of the benchmark
    executables (`ON` by default);
  * `DETRAY_BENCHMARKS_MULTITHREAD`: Boolean option making the benchmarks
    multithreaded (`OFF` by default);
  * `DETRAY_BENCHMARKS_REP`: String option with an integer for the repetitions
    that the benchmarks should run (`1` by default).

The following options configure how the build should set up the externals that
it needs:

- `DETRAY_SETUP_<XXX>`: Boolean to turn on/off the explicit "setup" of
  the externals (`ALGEBRA_PLUGINS`, `VECMEM`, `DFELIBS`, `MATPLOTPP`, `THRUST`,
  `GOOGLETEST` and `BENCHMARK`);
- `DETRAY_USE_SYSTEM_<XXX>`: Boolean configuring how to set up a given external
  * `ON`: The external is searched for "on the system" using
    [find_package](https://cmake.org/cmake/help/latest/command/find_package.html);
  * `OFF`: The package is set up for build as part of this project, using
    [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html).

### Benchmark Monitoring

A simple regression benchmark monitoring using `google/benchmark` runs on merge into master, the results are shown on https://acts-project.github.io/detray/.
