# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the project.
cmake_minimum_required( VERSION 3.11 )
project( detray VERSION 0.16 LANGUAGES CXX )

# Set up the used C++ standard(s).
set( CMAKE_CXX_STANDARD 17 CACHE STRING "The (host) C++ standard to use" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Disable (host) C++ extensions" )
set( CMAKE_CUDA_STANDARD 17 CACHE STRING "The (CUDA) C++ standard to use" )
set( CMAKE_CUDA_EXTENSIONS FALSE CACHE BOOL "Disable (CUDA) C++ extensions" )
set( CMAKE_SYCL_STANDARD 17 CACHE STRING "The (SYCL) C++ standard to use" )

# CMake include(s).
include( CMakeDependentOption )
include( GNUInstallDirs )

# Flags controlling the meta-build system.
option( DETRAY_USE_SYSTEM_LIBS "Use system libraries be default" FALSE )

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

# Include the Detray CMake code.
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )
include( detray-functions )

# Check if CUDA is available.
include( CheckLanguage )
check_language( CUDA )
set( DETRAY_BUILD_CUDA_DEFAULT FALSE )
if( CMAKE_CUDA_COMPILER )
   set( DETRAY_BUILD_CUDA_DEFAULT TRUE )
endif()

# Flags controlling which parts of Detray to build.
option( DETRAY_EIGEN_PLUGIN "Build Eigen math plugin" ON )
option( DETRAY_SMATRIX_PLUGIN "Build ROOT/SMatrix math plugin" OFF )
option( DETRAY_VC_PLUGIN "Build Vc based math plugin" ON )
option( DETRAY_IO_CSV "Build CSV IO module" ON )
option( DETRAY_IO_JSON "Build JSON IO module" ON)
option( DETRAY_SVG_DISPLAY "Build ActSVG display module" OFF)
option( DETRAY_BUILD_SYCL "Build the SYCL sources included in traccc" OFF )
option( DETRAY_BUILD_CUDA "Build the CUDA sources included in detray"
   ${DETRAY_BUILD_CUDA_DEFAULT} )
option( DETRAY_BUILD_TESTING "Build the (unit) tests of Detray"
   TRUE )
cmake_dependent_option( DETRAY_BENCHMARKS "Enable benchmark tests" TRUE
   "DETRAY_BUILD_TESTING" OFF )
cmake_dependent_option( DETRAY_BENCHMARKS_MULTITHREAD
   "Enable multithreaded benchmarks" OFF "DETRAY_BENCHMARKS" OFF )
set( DETRAY_BENCHMARKS_REP 1 CACHE STRING
   "Number of repetitions for the benchmarks" )

# Clean up.
unset( DETRAY_BUILD_CUDA_DEFAULT )

# Set up VecMem.
option( DETRAY_SETUP_VECMEM
   "Set up the VecMem target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_VECMEM
   "Pick up an existing installation of VecMem from the build environment"
   ${DETRAY_USE_SYSTEM_LIBS} )
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
   ${DETRAY_USE_SYSTEM_LIBS} )
if( DETRAY_SETUP_ALGEBRA_PLUGINS )
   if( DETRAY_USE_SYSTEM_ALGEBRA_PLUGINS )
      find_package( algebra-plugins REQUIRED )
   else()
      add_subdirectory( extern/algebra-plugins )
   endif()
endif()

# Set up ACTSVG for displaying
option( DETRAY_SETUP_ACTSVG
   "Set up the actsvg target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_ACTSVG
   "Pick up an existing installation of actsvg from the build environment"
   ${DETRAY_USE_SYSTEM_LIBS} )
if( DETRAY_SETUP_ACTSVG )
   if( DETRAY_USE_SYSTEM_ACTSVG )
      find_package( actsvg::core actsvg::meta REQUIRED )
   else()
      add_subdirectory( extern/actsvg )
   endif()
endif()

# Set up JSON for I/O
option( DETRAY_SETUP_NLOHMANN
   "Set up the nlohmann::json target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_NLOHMANN
   "Pick up an existing installation of nlohman::json from the build environment"
   ${DETRAY_USE_SYSTEM_LIBS} )
if( DETRAY_SETUP_NLOHMANN )
   if( DETRAY_USE_SYSTEM_NLOHMANN )
      find_package( nlohmann_json REQUIRED )
   else()
      add_subdirectory( extern/nlohmann_json )
   endif()
endif()

# Set up dfelibs.
option( DETRAY_SETUP_DFELIBS
   "Set up the dfelibs target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_DFELIBS
   "Pick up an existing installation of dfelibs from the build environment"
   ${DETRAY_USE_SYSTEM_LIBS} )
if( DETRAY_SETUP_DFELIBS )
   if( DETRAY_USE_SYSTEM_DFELIBS )
      find_package( dfelibs REQUIRED )
   else()
      add_subdirectory( extern/dfelibs )
   endif()
endif()

# Set up Thrust.
option( DETRAY_SETUP_THRUST
   "Set up the Thrust target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_THRUST
   "Pick up an existing installation of Thrust from the build environment"
   ${DETRAY_USE_SYSTEM_LIBS} )
if( DETRAY_SETUP_THRUST )
   if( DETRAY_USE_SYSTEM_THRUST )
      find_package( Thrust REQUIRED )
   else()
      add_subdirectory( extern/thrust )
   endif()
endif()
# Set up an IMPORTED library on top of the Thrust library/libraries. One that
# the Detray code could depend on publicly.
set( DETRAY_THRUST_OPTIONS "" CACHE STRING
   "Extra options for configuring how Thrust should be used" )
mark_as_advanced( DETRAY_THRUST_OPTIONS )
thrust_create_target( detray::Thrust ${DETRAY_THRUST_OPTIONS} )

# Set up GoogleTest.
option( DETRAY_SETUP_GOOGLETEST
   "Set up the GoogleTest target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_GOOGLETEST
   "Pick up an existing installation of GoogleTest from the build environment"
   ${DETRAY_USE_SYSTEM_LIBS} )
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
   ${DETRAY_USE_SYSTEM_LIBS} )
if( DETRAY_SETUP_BENCHMARK )
   if( DETRAY_USE_SYSTEM_BENCHMARK )
      find_package( benchmark REQUIRED )
   else()
      add_subdirectory( extern/benchmark )
   endif()
endif()

# Set up covfie.
option( DETRAY_SETUP_COVFIE
   "Set up the covfie target(s) explicitly" TRUE )
option( DETRAY_USE_SYSTEM_COVFIE
   "Pick up an existing installation of covfie from the build environment"
   ${DETRAY_USE_SYSTEM_LIBS} )
if( DETRAY_SETUP_COVFIE )
   if( DETRAY_USE_SYSTEM_COVFIE )
      find_package( covfie 0.5.0 REQUIRED )
   else()
      add_subdirectory( extern/covfie )
   endif()
endif()

# Set up all of the libraries of the project.
add_subdirectory( core )
add_subdirectory( plugins )
add_subdirectory( io )
add_subdirectory( utils )

# Set up the test(s).
include( CTest )
if( BUILD_TESTING AND DETRAY_BUILD_TESTING )
   add_subdirectory( tests )
endif()

# Set up the packaging of the project.
include( detray-packaging )
