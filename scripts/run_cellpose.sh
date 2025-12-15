#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <dir> <start> <end>"
    exit 1
fi

DIR=$1
RANGE_START=$2
RANGE_END=$3

for ((start=RANGE_START; start<RANGE_END; start+=1)); do
    end=$((start + 1))
    echo "Running with argument: $DIR and $start and $end"
    python run_cellpose.py "$DIR" "$start" "$end"
    echo "------------------------------------"
done
