# detray

Detray is part of the ACTS project (R&D line for parallelization), the ACTS project can be found: https://github.com/acts-project/acts.

This is a C++20 header only library for detector surface intersections using different algebra plugin libraries. It follows the navigation and propagation concept of ACTS, however, with an attempt to create
a geometry without polymorphic inheritance structure.


## Requirements and dependencies
#### OS & compilers:

- The C++ compiler must support C++20
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

A number of cmake preset configurations are provided and can be listed by:
```shell
cmake -S detray --list-presets
```
For a developer build, the ```dev-fp32``` and ```dev-fp64``` configurations are available (```fp```: floating point precision):
```shell
cmake -S detray -B detray-build --preset dev-fp32
```
The developer presets will fetch all dependencies, but not automatically trigger the build of additional detray components. For example, in order to trigger the build of the unit tests, the corresponding option needs to be specified:
```shell
cmake -S detray -B detray-build --preset dev-fp32 \
-DDETRAY_BUILD_UNITTESTS=ON
```
A full build, containing all components (e.g. tests and benchmarks), can be configured using the ```full-fp32``` and ```full-fp64``` presets.

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

Given a detector (and optionally also a grid and a material) json file, a number of validation test can be run from the command-line. For this, the library has to be built with the ```-DDETRAY_BUILD_CLI_TOOLS=ON``` flag. And example detector file can then be obtained using e.g.
```shell
detray-build/bin/detray_generate_toy_detector --write_material --write_grids
```
All of the validation tools presented in the following can also be run as part of a corresponding python script which takes the same arguments and will automatically create plots from the collected data. However, this requires Python 3, pandas, SciPy and NumPy, as well as Matplotlib to be available.

The detector geometry can be visualized in svg format with the following command:
```shell
detray-build/bin/detray_detector_display \
   --geometry_file  ./toy_detector/toy_detector_geometry.json [OPTION]...
```
The tool can also display single volumes or surfaces, as well as the navigation grids and material maps (the corresponding json files need to loaded in this case).

### Navigation Validation

In order to validate that the navigation works correctly in a given detector geometry, run the detector validation tool. It will first perform a consistency check on the detector, followed by a "ray scan" of the detector. The scan result will be compared to a full straight-line navigation run for every ray. After that, the navigation in a constant magnetic field of 2T is being tested in a similar fashion, using parameterized helix trajectories and a Newton-Raphson/Bisection algorithm to generate the truth intersections. For example:
```shell
detray-build/bin/detray_detector_validation \
    --geometry_file ./toy_detector/toy_detector_geometry.json \
    --grid_file ./toy_detector/toy_detector_surface_grids.json \
    --search_window 3 3 --n_tracks 100 --pT_range 0.5 100
```
In case of failures, this command will give a detailed debug output in the form of a log file, as well as an svg representation of the failed tracks. The grid file is optional, but will trigger the use of spacial grids as acceleration structures during the navigation run.

### Material Validation

This tool checks whether the navigator picks up the material correctly by comparing the material found during a ray scan with the material collected during navigation by a specialized actor:
```shell
detray-build/bin/detray_material_validation \
    --geometry_file toy_detector/toy_detector_geometry.json \
    --material_file toy_detector/toy_detector_homogeneous_material.json \
    --phi_steps 100 --eta_steps 100 --eta_range -4 4
```
Note: The correct material file must be loaded in addition to the geometry file!

## Continuous benchmark

Monitoring the propagation throughput with the toy geometry

<img src="https://gitlab.cern.ch/acts/detray-benchmark/-/raw/master/plots/array_data.png?ref_type=heads" />
