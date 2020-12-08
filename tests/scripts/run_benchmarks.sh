#!/bin/bash

echo "Running core.benchmarks ..."
./bin/core_masks

echo "Running eigen.benchmarks ..."
./bin/eigen_intersect_surfaces
./bin/eigen_intersect_all
