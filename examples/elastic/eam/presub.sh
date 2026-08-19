#!/bin/bash

# You can copy this to your submit script.

for dir in sample_*; do
    cd $dir
    echo "Running MD sample in $dir..."
    gpumd > log
    cd ..
done
