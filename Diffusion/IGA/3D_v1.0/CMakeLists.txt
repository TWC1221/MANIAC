cmake_minimum_required(VERSION 3.15)
project(neutron_code LANGUAGES Fortran C)

set(CMAKE_BUILD_TYPE "Debug")

# 1. HARD-CODED PATHS
set(PETSC_DIR "/usr/lib/petsc")
set(PETSC_INC "${PETSC_DIR}/include")
set(PETSC_LIB_DIR "${PETSC_DIR}/lib")

# 2. COMPILER FLAGS
# Use -I for modules and PETSc headers
set(CMAKE_Fortran_FLAGS "-Wall -ffree-line-length-none -I${PETSC_INC}")
set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -Og -fopenmp")
# -L tells the linker where to find the files named 'petsc', 'mpi', etc.
set(CMAKE_EXE_LINKER_FLAGS "-L${PETSC_LIB_DIR} -Wl,-rpath,${PETSC_LIB_DIR}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)

# 3. STATIC LIBRARY
add_library(src STATIC
    src/m_constants.F90
    src/m_types.F90
    src/m_asmg.F90
    src/m_quadrature.F90
    src/m_petsc.F90
    src/m_PCG.F90
    src/m_basis.F90
    src/m_material.F90
    src/m_outpVTK.F90
    src/m_boundaries.F90
)
set_target_properties(src PROPERTIES Fortran_PREPROCESS ON)

# 4. EXECUTABLE (Updated to 3D)
add_executable(3DIGA 3D_V1.0.F90)
set_target_properties(3DIGA PROPERTIES Fortran_PREPROCESS ON)

# 5. LINKING
target_link_libraries(3DIGA 
    src 
    petsc 
    mpifort 
    mpi 
    lapack 
    blas 
    m 
    dl
)