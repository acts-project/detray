#!/bin/bash
#
# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

echo "===> CI Benchmark running script for detray"

# Set the script path
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set workspace directory
if [ -z "${GITHUB_WORKSPACE}" ]; then
    WORKSPACE=${BASEDIR}/../../
else
    WORKSPACE=${GITHUB_WORKSPACE}
fi

git config --global --add safe.directory $(pwd)

# Generate a CSV file in data directory
export LASTCOMMIT=$(git log -n1 | head -n1 | cut -b 8-14)
export DETRAY_TEST_DATA_DIR=${WORKSPACE}/data/

touch ${WORKSPACE}/benchmark_${LASTCOMMIT}.csv

echo "===> Benchmark pipeline for commit ${LASTCOMMIT}"

# Run the benchmarks
for group in eigen array ; do
    echo "===> Running ${group}.benchmarks ..."
    ${WORKSPACE}/build/bin/detray_${group}_masks --benchmark_out=${group}_masks.csv --benchmark_out_format=csv
    ${WORKSPACE}/build/bin/detray_${group}_intersect_surfaces --benchmark_out=${group}_intersect_surfaces.csv --benchmark_out_format=csv
    ${WORKSPACE}/build/bin/detray_${group}_intersect_all --benchmark_out=${group}_intersect_all.csv --benchmark_out_format=csv

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

pip3 install matplotlib==3.7.1 numpy==1.24.2 pandas==1.5.3

echo "===> Download benchmark history ..."

# Configure Github user
git config --local user.email "action@github.com"
git config --local user.name "GitHub Action"

# Checkout to benchmark branch
git fetch origin
git checkout -b gh-pages origin/gh-pages

# Bring analysis script
git checkout ${LASTCOMMIT} -- ${WORKSPACE}/tests/scripts/analyze_benchmarks.py

# Record the benchmark result into benchmarks_history
BENCHMARK_HISTORY_FILE=${WORKSPACE}/archive/benchmarks/benchmarks_history.csv
cat ${WORKSPACE}/benchmark_${LASTCOMMIT}.csv >> ${BENCHMARK_HISTORY_FILE}

echo "===> Running benchmark analysis ..."

cd ${WORKSPACE}
python3 ${BASEDIR}/analyze_benchmarks.py ${BENCHMARK_HISTORY_FILE}

echo "===> Commiting the results ..."

cp ${WORKSPACE}/*.png  ${WORKSPACE}/figures/.

git add ${BENCHMARK_HISTORY_FILE}
git add ${WORKSPACE}/figures/*.png

git commit -m "updating benchmark data for commit ${LASTCOMMIT}" -a

if [ -z "${GITHUB_WORKSPACE}" ]; then
    echo "===> Exit"
else
    echo "===> Uploading the results"
    # CI stuck here...
    #git push origin gh-pages
    echo "===> Exit"
fi
