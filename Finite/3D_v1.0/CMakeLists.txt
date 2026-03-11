cmake_minimum_required(VERSION 3.1)
project(neutron_code LANGUAGES Fortran C)

set(CMAKE_BUILD_TYPE "Debug")

# 1. HARD-CODED PATHS (Matching your Ubuntu /usr/lib/petsc layout)
set(PETSC_DIR "/usr/lib/petsc")
set(PETSC_INC "${PETSC_DIR}/include")
set(PETSC_LIB_DIR "${PETSC_DIR}/lib")

# 2. COMPILER FLAGS
# We add the include path here and also tell the linker where the libraries live
set(CMAKE_Fortran_FLAGS "-Wall -ffree-line-length-none -I${PETSC_INC}")
set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -pg -Og -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "-L${PETSC_LIB_DIR} -Wl,-rpath,${PETSC_LIB_DIR}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# 3. STATIC LIBRARY (Modules)
add_library(src STATIC
    src/m_constants.F90
    src/m_types.F90
    src/m_quadrature.F90
    src/m_petsc.F90
    src/m_PCG.F90
    src/m_outpVTK.F90
    src/m_GMSH.F90
    src/m_finite_elements.F90
    src/m_material.F90
    src/m_boundaries.F90
    src/m_mesh_quality.F90
    src/m_timing.F90
)
set_target_properties(src PROPERTIES Fortran_PREPROCESS ON)

# 4. EXECUTABLE
add_executable(3Diffusion 3D_diffusion.F90)
set_target_properties(3Diffusion PROPERTIES Fortran_PREPROCESS ON)

# 5. THE LINKING FIX
# We link PETSc, and then explicitly link the MPI and Lapack/Blas libraries 
# that the DSO error was complaining about.
target_link_libraries(3Diffusion 
    src 
    petsc 
    mpifort 
    mpi 
    lapack 
    blas 
    m 
    dl
)




# cmake_minimum_required(VERSION 3.1)
# project(neutron_code LANGUAGES Fortran C)

# set(CMAKE_BUILD_TYPE "Debug")

# set(PETSC_INCLUDES "/usr/lib/petsc/include")
# set(PETSC_LIBS "/usr/lib/petsc/lib/libpetsc.so")

# set(CMAKE_Fortran_FLAGS "-Wall -ffree-line-length-none -I${PETSC_INCLUDES}")
# set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -pg -Og -fopenmp")

# set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)
# file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# add_library(src STATIC
#     src/m_constants.F90
#     src/m_types.F90
#     src/m_quadrature.F90
#     src/m_petsc.F90
#     src/m_PCG.F90
#     src/m_outpVTK.F90
#     src/m_GMSH.F90
#     src/m_finite_elements.F90
#     src/m_material.F90
#     src/m_boundaries.F90
#     src/m_mesh_quality.F90
# )

# set_target_properties(src PROPERTIES Fortran_PREPROCESS ON)

# add_executable(3Diffusion 3D_diffusion.F90)
# #add_executable(exe C5G7_MAKER.F90)
# set_target_properties(3Diffusion PROPERTIES Fortran_PREPROCESS ON)

# target_link_libraries(3Diffusion src ${PETSC_LIBS} lapack blas)