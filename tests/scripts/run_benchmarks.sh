#!/bin/bash

echo "===> CI Benchmark running script for detray"

export PWD_BUILD=$(pwd)
cd ${GITHUB_WORKSPACE}
export LASTCOMMIT=$(git log -n1 | head -n1 | cut -b 8-14)
echo "===> Benchmark pipeline for commit ${LASTCOMMIT}"
cd ${PWD_BUILD}

export DETRAY_TEST_DATA_DIR=${GITHUB_WORKSPACE}/tests/data

echo "===> Running eigen.benchmarks ..."
./bin/eigen_masks --benchmark_out=eigen_masks.csv --benchmark_out_format=csv
./bin/eigen_intersect_surfaces --benchmark_out=eigen_intersect_surfaces.csv --benchmark_out_format=csv
./bin/eigen_intersect_all --benchmark_out=eigen_intersect_all.csv --benchmark_out_format=csv

echo "===> Extracting benchmark results ..."
cat eigen_masks.csv | tail -n5  > eigen_masks_cropped.csv 
cat eigen_intersect_surfaces.csv | tail -f -n3 > eigen_intersect_surfaces_cropped.csv
cat eigen_intersect_all.csv | tail -f -n1 > eigen_intersect_all_cropped.csv
sed -i -e 's/"BM_/'$LASTCOMMIT',"eigen","BM_/g' eigen_masks_cropped.csv
sed -i -e 's/"BM_/'$LASTCOMMIT',"eigen","BM_/g' eigen_intersect_surfaces_cropped.csv
sed -i -e 's/"BM_/'$LASTCOMMIT',"eigen","BM_/g' eigen_intersect_all_cropped.csv
cat eigen_masks_cropped.csv > benchmark_${LASTCOMMIT}.csv
cat eigen_intersect_surfaces_cropped.csv >> benchmark_${LASTCOMMIT}.csv
cat eigen_intersect_all_cropped.csv >> benchmark_${LASTCOMMIT}.csv

cat benchmark_${LASTCOMMIT}.csv

echo "===> Install components for benchmark analysis ..."
pip3 install matplotlib numpy pandas

echo "===> Download benchmark history ..."

cd ${GITHUB_WORKSPACE} 
git fetch
git checkout -b gh-pages origin/gh-pages
cp archive/benchmarks/benchmarks_history.csv ${PWD_BUILD}/.
cd ${PWD_BUILD}
cat benchmark_${LASTCOMMIT}.csv >> benchmarks_history.csv

echo "===> Running benchmark analysis ..."
python3 -i analyze_benchmarks.py

echo "===> Prepare for uploading results ..."
cp benchmarks_history.csv ${GITHUB_WORKSPACE}/archive/benchmarks/benchmarks_history.csv
cp *.png  ${GITHUB_WORKSPACE}/figures/.
cd ${GITHUB_WORKSPACE}
git config --local user.email "action@github.com"
git config --local user.name "GitHub Action"
git add archive/benchmarks/benchmarks_history.csv
git add figures/*.png
git commit -m"updating benchmark data for commit ${LASTCOMMIT}" -a
