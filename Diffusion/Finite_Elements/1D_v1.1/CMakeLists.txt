cmake_minimum_required(VERSION 3.1) 

# Set our project name and Fortran as only language we compile 
# (C and C++ are also available) and should be added if we interop 
project(neutron_code VERSION 0.1 
    DESCRIPTION "PhD code"
    LANGUAGES Fortran C
    )

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# PETSc configuration
#set(PETSC_DIR "$ENV{PETSC_DIR}")
#set(PETSC_INCLUDES "${PETSC_DIR}/include")
#set(PETSC_LIBRARIES "${PETSC_DIR}/lib/libpetsc.so")

if(APPLE)
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -ffree-line-length-none ")
        set(CMAKE_Fortran_FLAGS_DEBUG   "-ffree-line-length-800 -O0 -fcheck=bounds -DDEBUG -fopenmp -llapack -lblas")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -DNDEBUG -fopenmp")
    endif()
else()
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall -ffree-line-length-none ")
	set(CMAKE_Fortran_FLAGS_DEBUG   "-fcheck=all -pedantic -ffree-line-length-800 -O0 -g -Og -fdec -DDEBUG -pg -fopenmp")
	set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -DNDEBUG -fopenmp")

    elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} ")
        set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -DDEBUG -pg -fopenmp")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -DNDEBUG -fopenmp")
    endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# Place Fortran module files in a known directory under the build tree so
# moving source files doesn't leave stale absolute paths to missing folders.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_mods)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

set(PETSC_DIR "$ENV{PETSC_DIR}")
set(PETSC_INCLUDES "${PETSC_DIR}/include")
set(PETSC_LIBRARIES "${PETSC_DIR}/lib/libpetsc.so")

message(STATUS "PETSc found in: ${PETSC_DIR}")
message(STATUS "PETSc include: ${PETSC_INCLUDES}")
message(STATUS "PETSC lib: ${PETSC_LIBRARIES}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${PETSC_DIR}/include")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -I${PETSC_DIR}/include")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -I${PETSC_DIR}/include")

# EXTERNAL PACKAGES
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(PkgConfig REQUIRED)

# Print PETSc details
# message(STATUS "Found PETSc: ${PETSC_DIR}")
# message(STATUS "PETSc Includes: ${PETSC_INCLUDES}")
# message(STATUS "PETSc Libraries: ${PETSC_LIBRARIES}")

# # Print SLEPc details
# message(STATUS "Found SLEPc: ${SLEPC_DIR}")
# message(STATUS "SLEPc Includes: ${SLEPC_INCLUDES}")
# message(STATUS "SLEPc Libraries: ${SLEPC_LIBRARIES}")

# Add link directories
add_library(src
src/m_constants.f90
src/m_quadrature.f90
src/m_materials.f90
src/m_petsc.F90
src/m_fem.F90
src/m_CSR.f90
src/m_boundaries.f90
src/m_assembly.F90
src/m_integrals.f90
src/m_PCG_solver.f90
src/m_output.f90
)

add_executable(exe src/FEM_1D.f90)

# Ensure targets search the generated module directory for .mod files
target_include_directories(src PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(exe PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Manually link PETSc
target_include_directories(exe PRIVATE ${PETSC_INCLUDES})
target_link_libraries(exe PRIVATE src ${PETSC_LIBRARIES} LAPACK::LAPACK BLAS::BLAS)
#target_include_directories(exe PRIVATE src)
#target_link_libraries(exe PRIVATE src LAPACK::LAPACK BLAS::BLAS)