#!/bin/bash

echo "Running core.benchmarks ..."
./tests/benchmarks/core/core_masks

echo "Running eigen.benchmarks ..."
./tests/benchmarks/eigen/eigen_intersect_surfaces
./tests/benchmarks/eigen/eigen_intersect_all
