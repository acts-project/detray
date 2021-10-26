# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the project.
cmake_minimum_required( VERSION 3.11 )
project( detray VERSION 0.1 LANGUAGES CXX )

# CMake include(s).
include( CMakeDependentOption )
include( GNUInstallDirs )

# Explicitly set the output directory for the binaries. Such that if this
# project is included by another project, the main project's configuration would
# win out.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" CACHE PATH
   "Directory for the built binaries" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built libraries" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built static libraries" )

# Set the language standards to use.
set( CMAKE_CXX_STANDARD 17 CACHE STRING "The (Host) C++ standard to use" )
set( CMAKE_CUDA_STANDARD 17 CACHE STRING "The (CUDA) C++ standard to use" )

option(DETRAY_ARRAY_PLUGIN "Build std::array math plugin" On)
option(DETRAY_VECMEM_ARRAY_PLUGIN "Build vecmem::array math plugin" On)
option(DETRAY_EIGEN_PLUGIN "Build Eigen math plugin" On)
option(DETRAY_SMATRIX_PLUGIN "Build ROOT/SMatrix math plugin" Off)
option(DETRAY_VC_ALGEBRA_PLUGIN "Build vc based simd algebra math plugin" ON)
option(DETRAY_IO_CSV "Build CSV IO module" On)
option(DETRAY_DISPLAY "Build matplot++ based display module" Off)
option(DETRAY_BUILD_CUDA "Build the CUDA sources included in detray" Off)
option(DETRAY_USE_VECMEM "Build with vecmem" On)
option( DETRAY_BUILD_TESTING "Build the (unit) tests of Detray"
   TRUE )
option( DETRAY_BENCHMARKS "Enable benchmark tests" TRUE )

# Set up VecMem.
option( DETRAY_SETUP_VECMEM
   "Set up the VecMem target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_VECMEM
   "Pick up an existing installation of VecMem from the build environment"
   FALSE )
if( DETRAY_SETUP_VECMEM )
   if( DETRAY_USE_SYSTEM_VECMEM )
      find_package( vecmem REQUIRED )
   else()
      add_subdirectory( extern/vecmem )
   endif()
endif()

# Set up Algebra Plugins.
option( DETRAY_SETUP_ALGEBRA_PLUGINS
   "Set up the Algebra Plugins target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_ALGEBRA_PLUGINS
   "Pick up an existing installation of Algebra Plugins from the build environment"
   FALSE )
if( DETRAY_SETUP_ALGEBRA_PLUGINS )
   if( DETRAY_USE_SYSTEM_ALGEBRA_PLUGINS )
      find_package( algebra-plugins REQUIRED )
   else()
      add_subdirectory( extern/algebra-plugins )
   endif()
endif()

# Set up dfelibs.
option( DETRAY_SETUP_DFELIBS
   "Set up the dfelibs target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_DFELIBS
   "Pick up an existing installation of dfelibs from the build environment"
   FALSE )
if( DETRAY_SETUP_DFELIBS )
   if( DETRAY_USE_SYSTEM_DFELIBS )
      find_package( dfelibs REQUIRED )
   else()
      add_subdirectory( extern/dfelibs )
   endif()
endif()

# Set up Matplot++.
option( DETRAY_SETUP_MATPLOTPP
   "Set up the Matplot++ target(s) explicitly" FALSE )
option( DETRAY_USE_SYSTEM_MATPLOTPP
   "Pick up an existing installation of Matplot++ from the build environment"
   FALSE )
if( DETRAY_SETUP_MATPLOTPP )
   if( DETRAY_USE_SYSTEM_MATPLOTPP )
      find_package( Matplot++ REQUIRED )
   else()
      add_subdirectory( extern/matplotplusplus )
   endif()
endif()

# Set up GoogleTest.
option( DETRAY_SETUP_GOOGLETEST
   "Set up the GoogleTest target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_GOOGLETEST
   "Pick up an existing installation of GoogleTest from the build environment"
   FALSE )
if( DETRAY_SETUP_GOOGLETEST )
   if( DETRAY_USE_SYSTEM_GOOGLETEST )
      find_package( GTest REQUIRED )
   else()
      add_subdirectory( extern/googletest )
   endif()
endif()

# Set up Google Benchmark.
option( DETRAY_SETUP_BENCHMARK
   "Set up the Google Benchmark target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_BENCHMARK
   "Pick up an existing installation of Google Benchmark from the build environment"
   FALSE )
if( DETRAY_SETUP_BENCHMARK )
   if( DETRAY_USE_SYSTEM_BENCHMARK )
      find_package( benchmark REQUIRED )
   else()
      add_subdirectory( extern/benchmark )
   endif()
endif()

cmake_dependent_option(DETRAY_BENCHMARKS_MULTITHREAD "Enable multithreaded benchmarks" Off "DETRAY_BENCHMARKS" Off)
if(DETRAY_BENCHMARKS_MULTITHREAD)
    message(STATUS "Using multithreaded benchmarks")
endif()
if(DETRAY_BENCHMARKS AND DETRAY_BENCHMARKS_REP)
    message(STATUS "Repetitions for benchmarks: " ${DETRAY_BENCHMARKS_REP})
endif()

# Set GPU-related options
if(DETRAY_BUILD_CUDA)
  enable_language(CUDA)
endif()

add_subdirectory(apps)
add_subdirectory(core)
add_subdirectory(plugins)
add_subdirectory(io)

# Set up the test(s).
include( CTest )
if( BUILD_TESTING AND DETRAY_BUILD_TESTING )
  add_subdirectory( tests )
endif()

if (DETRAY_DISPLAY)
    add_subdirectory(display)
endif()
