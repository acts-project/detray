#!/bin/bash

# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Number of threads
n_threads=1

# Number of tracks per thread
n_tracks_per_thread=1000

# Min and Max log10 rk tolerance to iterate
log10_min_rk_tol=-6
log10_max_rk_tol=2

# log10 rk tolerance for covariance transport
log10_rk_tol=-4

# Helix intersector tolerance [log10 in mm]
log10_helix_tol=-3
# Surface tolerance [log10 in mm]
log10_on_surface_tol=-3

# Monte-Carlo seed
mc_seed=0

# include multiple scattering
include_multiple_scattering=false

# Skip the first phase
skip_first_phase=false

# Skip the second phase
skip_second_phase=false

# Verbose level
verbose_level=1

while getopts "hd:n:t:p:q:c:i:s:m:f:r:v:" arg; do
    case $arg in
        h)
            echo ""
            echo "Mandatory arguments"
            echo "-d <Directory of detray_integration_test_jacobian_validation>"
            echo ""
            echo "Optonal arguments"
            echo "-n <Number of threads>"
            echo "-t <Number of tracks per thread>"
            echo "-p <log10(min_rk_error_tolerance_in_mm)>"
            echo "-q <log10(max_rk_error_tolerance_in_mm)>"
            echo "-c <log10(rk_error_tolerance_in_mm_for_covariance_transport)>"
            echo "-i <log10(intersection_tolerance_in_mm)>"
            echo "-s <Monte-Carlo seed>"
            echo "-m <Include myltiple scattering>"
            echo "-f <Skip the first phase>"
            echo "-r <Skip the second phase>"
            echo "-v <Verbose level>"
            echo ""
            exit 0
        ;;
        d)
            dir=$OPTARG
            echo "Directory of detray_integration_test_jacobian_validation: ${dir}"
        ;;
        n)
            n_threads=$OPTARG
            echo "Number of threads: ${n_threads}"
        ;;
        t)
            n_tracks_per_thread=$OPTARG
            echo "Number of tracks per thread: ${n_tracks_per_thread}"
        ;;
        p)
            log10_min_rk_tol=$OPTARG
            echo "log10(min_rk_error_tolerance_in_mm): ${log10_min_rk_tol}"
        ;;
        q)
            log10_max_rk_tol=$OPTARG
            echo "log10(max_rk_error_tolerance_in_mm): ${log10_max_rk_tol}"
        ;;
        c)
            log10_rk_tol=$OPTARG
            echo "log10(rk_error_tolerance_in_mm_for_covariance_transport): ${log10_rk_tol}"
        ;;
        i)
            log10_helix_tol=$OPTARG
            log10_on_surface_tol=$OPTARG
            echo "log10(intersection_tolerance_in_mm): ${log10_helix_tol}"
        ;;
        s)
            mc_seed=$OPTARG
            echo "Monte-Carlo seed: ${mc_seed}"
        ;;
        m)
            include_multiple_scattering=$OPTARG
            echo "Include multiple scattering in the second phase: ${include_multiple_scattering}"
        ;;
        f)
            skip_first_phase=$OPTARG
            echo "Skip the first phase: ${skip_first_phase}"
        ;;
        r)
            skip_second_phase=$OPTARG
            echo "Skip the second phase: ${skip_second_phase}"
        ;;
        v)
            verbose_level=$OPTARG
            echo "Set the verbose level: ${verbose_level}"
        ;;
    esac
done

echo ""

if [ -z "${dir}" ]; then
    echo "Option -d is missing"
    exit 1
fi

##########################
# RK tolerance iteration #
##########################

if [ "$skip_first_phase" = false ] ; then
    
    echo "Starting rk toleracne iteration..."
    
    # Remove the old directories
    for (( i=0; i < ${n_threads}; ++i ))
    do
        rm -rf ${PWD}/thread_${i}
    done
    rm -rf ${PWD}/merged
    
    for (( i=0; i < ${n_threads}; ++i ))
    do
        n_skips=`expr ${i} \* ${n_tracks_per_thread}`
        
        command_rk_tolerance="${dir}/detray_integration_test_jacobian_validation \
        --output-directory=thread_${i} \
        --rk-tolerance-iterate-mode=true \
        --n-tracks=${n_tracks_per_thread} \
        --n-skips=${n_skips} \
        --log10-min-rk-tolerance=${log10_min_rk_tol} \
        --log10-max-rk-tolerance=${log10_max_rk_tol} \
        --log10-helix-tolerance=${log10_helix_tol} \
        --log10-on-surface-tolerance=${log10_on_surface_tol} \
        --mc-seed=${mc_seed} \
        --include-multiple-scattering=${include_multiple_scattering} \
        --verbose-level=${verbose_level}"
        ${command_rk_tolerance} &
    done
    wait
    
    echo "Finished rk toleracne iteration"
    
fi

#####################################
# Jacobi validation & Cov transport #
#####################################

if [ "$skip_second_phase" = false ] ; then
    
    echo "Starting Jacobi validation & Cov transport..."
    
    for (( i=0; i < ${n_threads}; ++i ))
    do
        n_skips=`expr ${i} \* ${n_tracks_per_thread}`
        
        command_jacobi_validation="${dir}/detray_integration_test_jacobian_validation \
        --output-directory=thread_${i} \
        --rk-tolerance-iterate-mode=false \
        --n-tracks=${n_tracks_per_thread} \
        --n-skips=${n_skips} \
        --log10-rk-tolerance=${log10_rk_tol} \
        --log10-helix-tolerance=${log10_helix_tol} \
        --log10-on-surface-tolerance=${log10_on_surface_tol} \
        --mc-seed=${mc_seed} \
        --include-multiple-scattering=${include_multiple_scattering} \
        --verbose-level=${verbose_level}"
        ${command_jacobi_validation} &
    done
    wait
    
    echo "Finished Jacobi validation & Cov transport"
fi

###################
# Merge Csv files #
###################

echo "Starting merging Csv files..."

file_names=()

# Get the unique file names
echo ""
echo "/// Merged Csv file list ///"
for full_name in ./thread_0/*; do
    # Only take csv format files
    if [[ "$full_name" == *".csv" ]];then
        name=$(basename -- "$full_name")
        file_names+=(${name})
        echo $name
    fi
done
echo ""

output_dir=merged
mkdir -p ${output_dir}
# Merge the files
for name in "${file_names[@]}"
do
    arr=()
    for (( i=0; i < ${n_threads}; ++i ))
    do
        arr+=(thread_${i}/${name})
    done
    awk 'FNR==1 && NR!=1{next;}{print}' ${arr[@]} > ./${output_dir}/${name}
done

echo "Finished merging Csv files"

####################
# Run ROOT Scripts #
####################

cd ${output_dir}

if [ "$skip_first_phase" = false ] && [ "$skip_second_phase" = false ]; then
    
    # Run rk_tolerance_comparision.C
    root -q '../../../tests/validation/root/rk_tolerance_comparison.C+O('${log10_min_rk_tol}','${log10_max_rk_tol}')'
    
    # Run jacobian_comparison.C
    root -q -l ../../../tests/validation/root/jacobian_comparison.C+O
    
    # Run covariance_validation.C
    root -q -l ../../../tests/validation/root/covariance_validation.C+O
    
    elif [ "$skip_first_phase" = true ] && [ "$skip_second_phase" = false ]; then
    
    # Run covariance_validation.C
    root ../../../tests/validation/root/covariance_validation.C+O
    
    elif [ "$skip_first_phase" = false ] && [ "$skip_second_phase" = true ]; then
    
    # Run rk_tolerance_comparision.C
    root '../../../tests/validation/root/rk_tolerance_comparison.C+O('${log10_min_rk_tol}','${log10_max_rk_tol}')'
fi
