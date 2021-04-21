#!/bin/bash

echo "===> CI Benchmark running script for detray"

export PWD_BUILD=$(pwd)
cd ${GITHUB_WORKSPACE}
export LASTCOMMIT=$(git log -n1 | head -n1 | cut -b 8-14)
echo "===> Benchmark pipeline for commit ${LASTCOMMIT}"
cd ${PWD_BUILD}

export DETRAY_TEST_DATA_DIR=${GITHUB_WORKSPACE}/data/

touch benchmark_${LASTCOMMIT}.csv

for group in eigen array ; do
    echo "===> Running ${group}.benchmarks ..."
    ./bin/${group}_masks --benchmark_out=${group}_masks.csv --benchmark_out_format=csv
    ./bin/${group}_intersect_surfaces --benchmark_out=${group}_intersect_surfaces.csv --benchmark_out_format=csv
    ./bin/${group}_intersect_all --benchmark_out=${group}_intersect_all.csv --benchmark_out_format=csv

    echo "===> Extracting benchmark results ..."
    cat ${group}_masks.csv | tail -n5  > ${group}_masks_cropped.csv 
    cat ${group}_intersect_surfaces.csv | tail -f -n3 > ${group}_intersect_surfaces_cropped.csv
    cat ${group}_intersect_all.csv | tail -f -n1 > ${group}_intersect_all_cropped.csv
    sed -i -e 's/"BM_/'$LASTCOMMIT',"'$group'","BM_/g' ${group}_masks_cropped.csv
    sed -i -e 's/"BM_/'$LASTCOMMIT',"'$group'","BM_/g' ${group}_intersect_surfaces_cropped.csv
    sed -i -e 's/"BM_/'$LASTCOMMIT',"'$group'","BM_/g' ${group}_intersect_all_cropped.csv
    cat ${group}_masks_cropped.csv >> benchmark_${LASTCOMMIT}.csv
    cat ${group}_intersect_surfaces_cropped.csv >> benchmark_${LASTCOMMIT}.csv
    cat ${group}_intersect_all_cropped.csv >> benchmark_${LASTCOMMIT}.csv
done

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
