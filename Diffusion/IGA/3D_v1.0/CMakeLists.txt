cmake_minimum_required(VERSION 3.15)
project(neutron_code LANGUAGES Fortran C)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)

# ---- FIND MPI ----
find_package(MPI REQUIRED)

# ---- MANUALLY SET PETSC (Homebrew Fix) ----
set(PETSC_DIR "/opt/homebrew/opt/petsc")
set(PETSC_INCLUDES "${PETSC_DIR}/include")
# This finds the actual library file in the Homebrew directory
find_library(PETSC_LIB NAMES petsc PATHS "${PETSC_DIR}/lib" NO_DEFAULT_PATH)

if(NOT PETSC_LIB)
    message(FATAL_ERROR "PETSc library not found in ${PETSC_DIR}/lib")
endif()

# ---- COMPILER FLAGS ----
set(CMAKE_Fortran_FLAGS "-Wall -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -Og")

# ---- LIBRARY ----
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

target_include_directories(src PRIVATE ${PETSC_INCLUDES} ${MPI_Fortran_INCLUDE_DIRS})
set_target_properties(src PROPERTIES Fortran_PREPROCESS ON)

# ---- EXECUTABLE ----
add_executable(3DIGA 3D_V1.0.F90)
set_target_properties(3DIGA PROPERTIES Fortran_PREPROCESS ON)

# ADD THIS LINE: Tells the compiler where to find <petsc/finclude/...>
target_include_directories(3DIGA PRIVATE ${PETSC_INCLUDES} ${MPI_Fortran_INCLUDE_DIRS})

# ---- LINKING ----
# On macOS, we use the Accelerate framework for LAPACK/BLAS
find_library(LAPACK_LIB Accelerate)

target_link_libraries(3DIGA
    src
    MPI::MPI_Fortran
    ${PETSC_LIB}
    ${LAPACK_LIB}
)