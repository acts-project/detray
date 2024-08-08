# detray

Detray is part of the ACTS project (R&D line for parallelization), the ACTS project can be found: https://github.com/acts-project/acts.

This is a C++17 header only library for detector surface intersections using different algebra plugin libraries. It follows the navigation and propagation concept of ACTS, however, with an attempt to create
a geometry without polymorphic inheritance structure.


## Requirements and dependencies
#### OS & compilers:

- The C++ compiler must support C++17
- The CUDA Toolkit version must be greater than major version 11

#### Dependency
- CMake

## Getting started

The respository is meant to be possible to build "out of the box", with standard
CMake build procedures.

```shell
git clone https://github.com/acts-project/detray.git
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -S detray -B detray-build
cmake --build detray-build
```

For tests and benchmarks with the inhomogeneous magnetic field, the ODD field in covfie format should be downloaded and the environment variable should be set.
```shell
cd detray/data
bash detray_data_get_files.sh
export DETRAY_BFIELD_FILE="${PWD}/odd-bfield_v0_9_0.cvf"
```

#### Build options

| Option | Description | Default |
| --- | --- | --- |
| DETRAY_BUILD_CUDA  | Build the CUDA sources included in detray | ON (if available) |
| DETRAY_BUILD_SYCL  | Build the SYCL sources included in detray | OFF |
| DETRAY_BUILD_TESTING  | Build the (unit) tests of detray | ON |
| DETRAY_BUILD_TUTORIALS  | Build the examples of detray | ON |
| DETRAY_CUSTOM_SCALARTYPE | Floating point precision | double |
| DETRAY_EIGEN_PLUGIN | Build Eigen math plugin | ON |
| DETRAY_SMATRIX_PLUGIN | Build ROOT/SMatrix math plugin | OFF |
| DETRAY_VC_PLUGIN | Build Vc based math plugin | ON |
| DETRAY_SVG_DISPLAY | Build ActSVG display module | OFF |

## Continuous benchmark

Monitoring the propagation speed with the toy geometry

<img src="https://gitlab.cern.ch/beyeo/detray-benchmark/-/raw/master/plots/array_data.png?ref_type=heads" width="500" height="500" /> 
<img src="https://gitlab.cern.ch/beyeo/detray-benchmark/-/raw/master/plots/eigen_data.png?ref_type=heads" width="500" height="500" />
