style=$1
export OPENBLAS_NUM_THREADS=1

if [ "$style" == "debug" ]; then
   echo "makeing clean debug"
   cmake -B build -DCMAKE_BUILD_TYPE='Debug'
   cmake --build build
elif [ "$style" == "release" ]; then
   echo "makeing clean release"
   cmake -B build -DCMAKE_BUILD_TYPE='Release'
   cmake --build build
else
   echo "makeing clean debug"
   cmake -B build -DCMAKE_BUILD_TYPE='Debug'
   cmake --build build
fi