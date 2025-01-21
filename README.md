# detray

[![Test Status](https://github.com/acts-project/detray/actions/workflows/builds.yml/badge.svg?branch=main)](https://github.com/acts-project/detray/actions/workflows/builds.yml)
[![Lint Status](https://github.com/acts-project/detray/actions/workflows/checks.yml/badge.svg?branch=main)](https://github.com/acts-project/detray/actions/workflows/checks.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=acts-project_detray&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=acts-project_detray)

Detray is part of the [ACTS project](https://github.com/acts-project/acts) (R&D line for parallelization).

This is a C++20 header-only library providing a GPU-friendly tracking detector description using different [linear algebra](https://github.com/acts-project/algebra-plugins) libraries. It follows the navigation and propagation concept of ACTS, implementing a geometry using a flat memory layout and no abstract interfaces (virtual functions). A detray detector can therefore be constructed on the host and copied to an accelerator device in a straight-forward way.

With the geometry description comes a fully featured, GPU-ready track state propagation implementation in inhomogeneous magnetic fields (vector field description using [covfie](https://github.com/acts-project/covfie)), with track parameter covariance transport including material interactions.

## Requirements and Dependencies
#### OS & Compilers:

- The C++ compiler must support C++20
- The CUDA Toolkit version must be greater than major version 11

#### Dependencies:
- CMake (version >= 3.14, version >= 3.18 for CUDA)


## Getting started

The respository should build "out of the box", with standard
CMake build procedures.

```shell
git clone https://github.com/acts-project/detray.git
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -S detray -B detray-build
cmake --build detray-build
```

For unit and integration tests using the *Open Data Detector* (ODD) solenoid field, a magnetic field map file in covfie format needs to be downloaded and the corresponding environment variable should be set to:
```shell
cd detray/data
bash detray_data_get_files.sh
export DETRAY_BFIELD_FILE="${PWD}/odd-bfield_v0_9_0.cvf"
```

#### Build options

A number of cmake preset configurations are provided and can be listed by:
```shell
cmake -S detray --list-presets
```
For a developer build, the `dev-fp32` and `dev-fp64` configurations are available (`fp`: floating point precision):
```shell
cmake -S detray -B detray-build --preset dev-fp32
```
The developer presets will fetch all dependencies, but not automatically trigger the build of additional detray components. For example, in order to trigger the build of the unit tests, the corresponding option needs to be specified:
```shell
cmake -S detray -B detray-build --preset dev-fp32 \
-DDETRAY_BUILD_UNITTESTS=ON
```
A full build, containing all components (e.g. tests and benchmarks), can be configured using the `full-fp32` and `full-fp64` presets.

The following cmake options are available and can also be specified explicitly for any preset:

| Option | Description | Default |
| --- | --- | --- |
| DETRAY_BUILD_CUDA  | Build the CUDA sources included in detray | ON (if available) |
| DETRAY_BUILD_SYCL  | Build the SYCL sources included in detray | OFF |
| DETRAY_BUILD_TEST_UTILS  | Build the detray test utilities library (contains e.g. test detectors) | OFF |
| DETRAY_BUILD_UNITTESTS  | Build the detray unit tests | OFF |
| DETRAY_BUILD_INTEGRATIONTESTS  | Build the detray integration tests | OFF |
| DETRAY_BUILD_ALL_TESTS  | Build the detray unit and integration tests | OFF |
| DETRAY_BUILD_BENCHMARKS  | Build the detray benchmarks | OFF |
| DETRAY_BUILD_CLI_TOOLS  | Build the detray command line tools | OFF |
| DETRAY_BUILD_TUTORIALS  | Build the examples of detray | OFF |
| DETRAY_CUSTOM_SCALARTYPE | Floating point precision | double |
| DETRAY_EIGEN_PLUGIN | Build Eigen math plugin | OFF |
| DETRAY_SMATRIX_PLUGIN | Build ROOT/SMatrix math plugin | OFF |
| DETRAY_VC_AOS_PLUGIN | Build Vc based AoS math plugin | OFF |
| DETRAY_VC_SOA_PLUGIN | Build Vc based SoA math plugin (currently only supports the ray-surface intersectors) | OFF |
| DETRAY_SVG_DISPLAY | Build ActSVG display module | OFF |


## Detector Validation

Given a detray detector (and optionally also a grid and a material) json file, a number of validation test can be run from the command-line. For this, the library has to be built with the `-DDETRAY_BUILD_CLI_TOOLS=ON` option enabled. An example detector file can then be obtained using e.g.
```shell
detray-build/bin/detray_generate_toy_detector --write_material --write_grids
```
All of the validation tools presented in the following can also be run as part of a corresponding [python script](https://github.com/acts-project/detray/tree/main/tests/tools/python) which takes the same arguments and will automatically create plots from the collected data. However, this requires Python 3, pandas, SciPy and NumPy, as well as Matplotlib to be available.

The detector geometry can be visualized in SVG format with the following command:
```shell
detray-build/bin/detray_detector_display \
   --geometry_file  ./toy_detector/toy_detector_geometry.json
```
The tool can also display single volumes or surfaces, as well as the navigation grids and material maps (the corresponding json files need to loaded in this case). For an overview of all available options for the command-line tools add `--help`.

### Navigation Validation

In order to validate that the navigation works correctly in a given detector geometry, run the detector validation tool. It will first perform a consistency check on the detector, followed by a "ray scan" of the detector. The scan result will be compared to a full straight-line navigation run for every ray. After that, the navigation in a constant magnetic field of 2T is being tested in a similar fashion, using parameterized helix trajectories and a Newton-Raphson/Bisection algorithm to generate the truth intersections. For example:
```shell
detray-build/bin/detray_detector_validation \
    --geometry_file ./toy_detector/toy_detector_geometry.json \
    --grid_file ./toy_detector/toy_detector_surface_grids.json \
    --search_window 3 3 --n_tracks 100 --pT_range 0.5 100
```
In case of failures, this command will give a detailed debug output in the form of a log file, as well as an SVG representation of the failed tracks. The grid file is optional, but will trigger the use of spacial grids as acceleration structures during the navigation run.

### Material Validation

This tool checks whether the navigator picks up the material correctly by comparing the material found during a ray scan with the material collected during navigation by a specialized actor:
```shell
detray-build/bin/detray_material_validation \
    --geometry_file toy_detector/toy_detector_geometry.json \
    --material_file toy_detector/toy_detector_homogeneous_material.json \
    --phi_steps 100 --eta_steps 100 --eta_range -4 4
```
Note: The correct material file must be loaded in addition to the geometry file!


## Benchmarks

A number of benchmarks exist, which are based on the google benchmark library, and can be run from command-line. For this, the `-DDETRAY_BUILD_BENCHMARKS=ON` and `-DDETRAY_BUILD_CLI_TOOLS=ON` flags need to be specified. Then pass the detray detector file(s) and additional options to the benchmark tools for the different hardware backends:
```shell
detray-build/bin/detray_propagation_benchmark_<backend>_<algebra> \
    --geometry_file ./toy_detector/toy_detector_geometry.json \
    --grid_file ./toy_detector/toy_detector_surface_grids.json \
    --material_file ./toy_detector/toy_detector_homogeneous_material.json \
    --sort_tracks --randomize_charge --eta_range -3 3 -pT_range 1 100
```
For every algebra-plugin that was built, a corresponding benchmark executable will be present. The CPU-backend benchmark is built by default and the CUDA-backend benchmark will be available if detray was built with CUDA enabled (`-DDETRAY_BUILD_CUDA=ON`).

### Continuous benchmark

Monitoring the propagation throughput with the toy geometry per commit:

<img src="https://gitlab.cern.ch/acts/detray-benchmark/-/raw/master/plots/array_data.png?ref_type=heads" />
