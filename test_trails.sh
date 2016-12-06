#!/bin/bash

set -o pipefail
set -o nounset
set -o errexit

MAT_SIZE=(1001 10001 100001 1000001 10000001 20000001 100000001)
THREAD_NUM=(1 2 3 4)

for MAT in "${MAT_SIZE[@]}"; do
    for THREADS in "${THREAD_NUM[@]}"; do
        gsed -Ei "s/#define NSUB .+/#define NSUB ${MAT}/" fem1d.c
        gsed -Ei "s/#define NUM_THREADS .+/#define NUM_THREADS ${THREADS}/" fem1d.c
        make clean
        make
        echo "Doing $THREADS $MAT"
        ./fem1d >> trails.csv
        ./fem1d >> trails.csv
        ./fem1d >> trails.csv
    done
done
