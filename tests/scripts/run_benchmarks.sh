#!/bin/bash

export lastcommit =`git log -n1 | head -n1 | cut -b 8-14`
echo "Benchmark nalysis of commit ${lastcommit}"

echo "Running core.benchmarks ..."
./bin/core_masks --benchmark_out=core_masks.csv --benchmark_out_format=csv

echo "Running eigen.benchmarks ..."
./bin/eigen_intersect_surfaces --benchmark_out=eigen_intersect_surfaces.csv --benchmark_out_format=csv
./bin/eigen_intersect_all  --benchmark_out=eigen_intersect_all.csv --benchmark_out_format=csv

echo "Install components for benchmark analysis ..."
pip3 install matplotlib numpy pandas


