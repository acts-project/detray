#!/bin/bash

lscpu

# Setup
module load gcc/13.2 cuda/12.6 cmake/3.30 clang-format/10 boost/1.87

# Build detray
#cd
#git clone https://github.com/acts-project/detray.git detray_sapphirerapids
#mkdir detray_sapphirerapids/build
#cd detray_sapphirerapids/build
#make clean

#cmake -S .. -B . --preset dev-fp32 -DCMAKE_BUILD_TYPE="Release" -DCMAKE_CXX_FLAGS="-march=native" \
#-DDETRAY_BUILD_TESTING=OFF -DDETRAY_BUILD_TUTORIALS=OFF -DDETRAY_BUILD_BENCHMARKS=ON -DDETRAY_BUILD_CLI_TOOLS=ON \
#-DDETRAY_EIGEN_PLUGIN=ON -DDETRAY_VC_AOS_PLUGIN=ON -DDETRAY_VC_SOA_PLUGIN=OFF \
#-DDETRAY_BUILD_CUDA=ON

make -j 16

# Prepwork

# Absolute path of detray git root directory
root_dir="$(git rev-parse --show-toplevel)"

# Assume location of installation
bin_dir="$root_dir/build/bin"

# Name of the hardware backend that the benchmarks run on
hw_backend="AMD EPYC 7413"
#hw_backend="RTX5000"
hw_backend_file="${hw_backend// /_}"

# How many times to repeat the benchmarks
repetitions=10

# Toy detector json files
toy_detector_geometry_file="$root_dir/build/toy_detector/toy_detector_geometry.json"
toy_detector_grid_file="$root_dir/build/toy_detector/toy_detector_surface_grids.json"
toy_detector_material_file="$root_dir/build/toy_detector/toy_detector_homogeneous_material.json"

# Open Data Detector (ODD) files
odd_geometry_file="/mnt/ssd1/jonierma/detray_new/data/odd-detray_geometry_detray.json"
odd_grid_file="/mnt/ssd1/jonierma/detray_new/data/odd-detray_surface_grids_detray.json"
odd_material_file="/mnt/ssd1/jonierma/detray_new/data/odd-detray_material_detray.json"

# Open Data Detector (ITk) files
itk_geometry_file="/mnt/ssd1/jonierma/detray_new/data/itk_geometry_detray.json"
itk_grid_file="/mnt/ssd1/jonierma/detray_new/data/itk_surface_grids_detray.json"
itk_material_file="/mnt/ssd1/jonierma/detray_new/data/itk_material_detray.json"

backend="cpu"
#backend="cuda"

n_tracks=10000

# Run all benchmarks
for algebra in array eigen vc_aos
do
    for pT in 1 10 100
    do
        for bench in benchmark scaling
        do

        # Run ITk benchmarks
        $bin_dir/detray_propagation_"$bench"_"$backend"_$algebra \
            --geometry_file $itk_geometry_file             \
            --n_tracks $n_tracks                           \
            --overstep_tolerance -300                      \
            --sort_tracks                                  \
            --randomize_charge                             \
            --eta_range -4 4                               \
            --pT_range $pT                                 \
            --bknd_name "$hw_backend"                      \
            --benchmark_repetitions=$repetitions           \
            --benchmark_display_aggregates_only=true       \
            --benchmark_out_format=json                    \
        --benchmark_filter=BM_PROPAGATION_ITk_250000_TRACKS*   \
            --benchmark_out=itk_"$bench"_"$backend"_"$hw_backend_file"_"$algebra"_pT"$pT"GeV.json

        $bin_dir/detray_propagation_"$bench"_"$backend"_$algebra \
            --geometry_file $itk_geometry_file             \
            --grid_file $itk_grid_file                     \
            --n_tracks $n_tracks                           \
            --overstep_tolerance -300                      \
            --sort_tracks                                  \
            --randomize_charge                             \
            --eta_range -4 4                               \
            --pT_range $pT                                 \
            --search_window 0 0                            \
            --bknd_name "$hw_backend"                      \
            --benchmark_repetitions=$repetitions           \
            --benchmark_display_aggregates_only=true       \
            --benchmark_out_format=json                    \
            --benchmark_out=itk_"$bench"_"$backend"_"$hw_backend_file"_"$algebra"_grids_pT"$pT"GeV.json

        $bin_dir/detray_propagation_"$bench"_"$backend"_$algebra \
            --geometry_file $itk_geometry_file             \
            --material_file $itk_material_file             \
            --n_tracks $n_tracks                           \
            --overstep_tolerance -300                      \
            --sort_tracks                                  \
            --randomize_charge                             \
            --covariance_transport                         \
            --eta_range -4 4                               \
            --pT_range $pT                                 \
            --bknd_name "$hw_backend"                      \
            --benchmark_repetitions=$repetitions           \
            --benchmark_display_aggregates_only=true       \
            --benchmark_out_format=json                    \
            --benchmark_out=itk_"$bench"_"$backend"_"$hw_backend_file"_"$algebra"_cov_pT"$pT"GeV.json

        $bin_dir/detray_propagation_"$bench"_"$backend"_$algebra \
            --geometry_file $itk_geometry_file             \
            --material_file $itk_material_file             \
            --n_tracks $n_tracks                           \
            --overstep_tolerance -300                      \
            --grid_file $itk_grid_file                     \
            --sort_tracks                                  \
            --randomize_charge                             \
            --covariance_transport                         \
            --eta_range -4 4                               \
            --pT_range $pT                                 \
            --search_window 0 0                            \
            --bknd_name "$hw_backend"                      \
            --benchmark_repetitions=$repetitions           \
            --benchmark_display_aggregates_only=true       \
            --benchmark_out_format=json                    \
            --benchmark_out=itk_"$bench"_"$backend"_"$hw_backend_file"_"$algebra"_cov_grids_pT"$pT"GeV.json
        done
    done
done
