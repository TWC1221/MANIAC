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
#     src/m_GMSH.F90
#     src/m_finite_elements.F90
#     src/m_material.F90
#     src/m_outpVTK.F90
#     src/m_boundaries.F90
# )

# set_target_properties(src PROPERTIES Fortran_PREPROCESS ON)

# add_executable(2Diffusion 2D_V1.1.F90)
# #add_executable(exe VTK_Extractor_V2.0.F90)
# set_target_properties(2Diffusion PROPERTIES Fortran_PREPROCESS ON)

# target_link_libraries(2Diffusion src ${PETSC_LIBS} lapack blas)

cmake_minimum_required(VERSION 3.5)
project(neutron_code LANGUAGES Fortran C)

set(CMAKE_BUILD_TYPE "Debug")

# 1. EXTERNAL PACKAGES
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# 2. COMPILER FLAGS
set(CMAKE_Fortran_FLAGS "-Wall -ffree-line-length-none")
if(APPLE)
    # -pg is often unsupported on macOS and can cause linker errors
    set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -Og -fopenmp")
else()
    set(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all -pg -Og -fopenmp")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)

# 3. STATIC LIBRARY (Modules)
add_library(src STATIC
    src/m_constants.F90
    src/m_types.F90
    src/m_quadrature.F90
    src/m_asmg.F90
    src/m_basis.F90
    src/m_material.F90
    src/m_precompute.F90
    src/m_outpVTK_dfem.F90
    src/m_sweep_order.F90
    src/m_transport.F90
)
set_target_properties(src PROPERTIES Fortran_PREPROCESS ON)

# 4. EXECUTABLE
add_executable(2DIGA 2D_V1.1.F90)
set_target_properties(2DIGA PROPERTIES Fortran_PREPROCESS ON)

# Ensure the executable finds the generated .mod files
target_include_directories(src PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(2DIGA PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# 5. THE LINKING FIX
target_link_libraries(2DIGA
    src
    LAPACK::LAPACK
    BLAS::BLAS
)