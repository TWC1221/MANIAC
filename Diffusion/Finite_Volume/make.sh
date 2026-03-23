style=$1
export OPENBLAS_NUM_THREADS=1

if [ "$style" == "debug" ]; then
   echo "-- Makeing clean debug..."
   echo "-- Removing old build directory..."
   rm -rf build
   cmake -B build -DCMAKE_BUILD_TYPE='Debug' -DCMAKE_Fortran_COMPILER=gfortran
   cmake --build build
elif [ "$style" == "release" ]; then
   echo "-- Makeing clean release..."
   echo "-- Removing old build directory..."
   rm -rf build 
   cmake -B build -DCMAKE_BUILD_TYPE='Release' -DCMAKE_Fortran_COMPILER=gfortran
   cmake --build build
else
   echo "-- Makeing clean debug..."
   echo "-- Removing old build directory..."
   rm -rf build
   cmake -B build -DCMAKE_BUILD_TYPE='Debug' -DCMAKE_Fortran_COMPILER=gfortran
   cmake --build build
fi

# style=$1           # e.g., debug or release
# compiler=$2        # e.g., ifx or gfortran

# export OPENBLAS_NUM_THREADS=1

# # Set default compiler if none provided
# if [ -z "$compiler" ]; then
#    compiler=gfortran
# fi

# if [ "$style" == "debug" ]; then
#    echo "making clean debug with $compiler"
#    rm -rf build
#    cmake -B build -DCMAKE_BUILD_TYPE='Debug' -DCMAKE_Fortran_COMPILER="$compiler"
#    cmake --build build
# elif [ "$style" == "release" ]; then
#    echo "making clean release with $compiler"
#    rm -rf build
#    cmake -B build -DCMAKE_BUILD_TYPE='Release' -DCMAKE_Fortran_COMPILER="$compiler"
#    cmake --build build
# else
#    echo "making clean debug with $compiler"
#    rm -rf build
#    cmake -B build -DCMAKE_BUILD_TYPE='Debug' -DCMAKE_Fortran_COMPILER="$compiler"
#    cmake --build build
# fi
