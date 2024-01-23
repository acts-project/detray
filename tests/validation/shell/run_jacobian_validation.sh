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

# Helix intersector tolerance [log10 in mm]
log10_helix_tol=-3
# Surface tolerance [log10 in mm]
log10_on_surface_tol=-3

while getopts "hd:n:t:p:q:i:" arg; do
    case $arg in
        h)
            echo ""
            echo "Mandatory arguments"
            echo "-d <Directory of detray_test_jacobian_validation>"
            echo ""
            echo "Optonal arguments"
            echo "-n <Number of threads>"
            echo "-t <Number of tracks per thread>"
            echo "-p <log10(min_rk_error_tolerance_in_mm)>"
            echo "-q <log10(max_rk_error_tolerance_in_mm)>"
            echo "-i <log10(intersection_tolerance_in_mm)>"
            echo ""
            exit 0
        ;;
        d)
            dir=$OPTARG
            echo "Directory of detray_test_jacobian_validation: ${dir}"
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
        i)
            log10_helix_tol=$OPTARG
            log10_on_surface_tol=$OPTARG
            echo "log10(intersection_tolerance_in_mm): ${log10_helix_tol}"
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

echo "Starting rk toleracne iteration..."

for (( i=0; i < ${n_threads}; ++i ))
do
    n_skips=`expr ${i} \* ${n_tracks_per_thread}`
    
    command_rk_tolerance="${dir}/detray_test_jacobian_validation \
    --output-directory=${i} \
    --rk-tolerance-iterate-mode=true \
    --n-tracks=${n_tracks_per_thread} \
    --n-skips=${n_skips} \
    --log10-min-rk-tolerance=${log10_min_rk_tol} \
    --log10-max-rk-tolerance=${log10_max_rk_tol} \
    --log10-helix-tolerance=${log10_helix_tol} \
    --log10-on-surface-tolerance=${log10_on_surface_tol}"
    
    ${command_rk_tolerance} &
done
wait

echo "Finished rk toleracne iteration"

#####################################
# Jacobi validation & Cov transport #
#####################################

echo "Starting Jacobi validation & Cov transport..."

for (( i=0; i < ${n_threads}; ++i ))
do
    n_skips=`expr ${i} \* ${n_tracks_per_thread}`
    
    command_jacobi_validation="${dir}/detray_test_jacobian_validation \
    --output-directory=${i} \
    --rk-tolerance-iterate-mode=false \
    --n-tracks=${n_tracks_per_thread} \
    --n-skips=${n_skips} \
    --log10-rk-tolerance=${log10_min_rk_tol} \
    --log10-helix-tolerance=${log10_helix_tol} \
    --log10-on-surface-tolerance=${log10_on_surface_tol}"
    
    ${command_jacobi_validation} &
done
wait

echo "Finished Jacobi validation & Cov transport"

###################
# Merge Csv files #
###################

echo "Starting merging Csv files..."

file_names=()

# Get the unique file names
for full_name in ./0/*; do
    name=$(basename -- "$full_name")
    file_names+=(${name})
done

output_dir=merged
mkdir -p ${output_dir}
# Merge the files
for name in "${file_names[@]}"
do
    arr=()
    for (( i=0; i < ${n_threads}; ++i ))
    do
        arr+=(${i}/${name})
    done
    awk 'FNR==1 && NR!=1{next;}{print}' ${arr[@]} > ./${output_dir}/${name}
done

echo "Finished merging Csv files"

####################
# Run ROOT Scripts #
####################

cd ${output_dir}

# Run rk_tolerance_comparision.C
root -q '../../../../tests/validation/root/rk_tolerance_comparison.C+O('${log10_min_rk_tol}','${log10_max_rk_tol}')'

# Run jacobian_comparison.C
root -q -l ../../../../tests/validation/root/jacobian_comparison.C+O

# Run covariance_validation.C
root -q -l ../../../../tests/validation/root/covariance_validation.C+O