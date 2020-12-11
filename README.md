# detray

Detray is part of the ACTS project (R&D line for parallelization), the ACTS project can be found: https://github.com/acts-project/acts.

This is a C++17 header only library for detector surface intersections using different algebra plugin libraries. It follows the navigation and propagation concept of ACTS, however, with an attempt to create
a geometry without polymorphic inheritance structure.

## Repository structure

The code for defining detector surfaces is within the `core` directory, the algebra plugins are in the `plugin` directory, and tests are all within the `test` library.

Testing and benchmarking is done with `googletest` and `google/benchmark`.


## Getting started

Clone the repository and initialize the submodules (googletest, benchmark).

```shell
git clone git@github.com:asalzburger/detray.git
cd detray
git submodule update --init
```

Running `CMake` (minimal version)
```shell
cmake -S . -B <build_directory>
```

Build
```shell
cmake --build <build_directory>
```

### Plugins

The following algebra plugins are avaiable and can be switched on/off 

| *Plugin* | *CMake Flag* | *Default value* | *Dependency requirement* |
| ---------|--------------|-----------------|--------------------------|
| array | DETRAY_ARRAY_PLUGIN | On | - |
| eigen | DETRAY_EIGEN_PLUGIN | Off | Eigen, http://eigen.tuxfamily.org/ |




### Benchmark Monitoring

A simple regression benchmark monitoring using `google/benchmark` runs on merge into master, the results are shown on https://asalzburger.github.io/detray/.
