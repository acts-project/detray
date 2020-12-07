#!/bin/bash

./tests/benchmarks/eigen/eigen_intersect_surfaces --benchmark_format=json | tee eigen_intersect_surfaces.json
./tests/benchmarks/eigen/eigen_intersect_all --benchmark_format=json | tee eigen_intersect_all.json
