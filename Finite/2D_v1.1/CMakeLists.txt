cmake_minimum_required(VERSION 3.1)
project(neutron_code LANGUAGES Fortran C)

set(CMAKE_BUILD_TYPE "Debug")

set(PETSC_INCLUDES "/usr/lib/petsc/include")
set(PETSC_LIBS "/usr/lib/petsc/lib/libpetsc.so")

set(CMAKE_Fortran_FLAGS "-Wall -ffree-line-length-none -I${PETSC_INCLUDES}")
set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -pg -Og -fopenmp")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(src STATIC
    src/m_constants.F90
    src/m_types.F90
    src/m_quadrature.F90
    src/m_petsc.F90
    src/m_PCG.F90
    src/m_GMSH.F90
    src/m_finite_elements.F90
    src/m_material.F90
    src/m_outpVTK.F90
    src/m_boundaries.F90
)

set_target_properties(src PROPERTIES Fortran_PREPROCESS ON)

add_executable(exe 2D_V1.1.F90)
#add_executable(exe VTK_Extractor_V2.0.F90)
set_target_properties(exe PROPERTIES Fortran_PREPROCESS ON)

target_link_libraries(exe src ${PETSC_LIBS} lapack blas)