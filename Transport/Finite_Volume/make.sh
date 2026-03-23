#!/bin/sh
style=$1
export OPENBLAS_NUM_THREADS=1
export PETSC_DIR=/usr/lib/petsc

if [ "$style" = "debug" ]; then
echo "Debug Version"
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
elif [ "$style" = "release" ]; then
echo "Release Version"
rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release  ..
else
echo "Debug Version"
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
fi
echo "Making..."
make -j

