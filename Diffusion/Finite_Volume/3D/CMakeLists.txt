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
src/m_CSR_types.f90
src/m_assembly.f90
src/m_subcritical.f90
src/m_BC_types.f90
src/m_current.f90
src/m_iter.F90
src/m_mesh_types.f90
src/m_outpVTK.f90
src/m_petsc.F90
src/m_PCG_solver.f90
src/m_wall_current.f90
src/m_kinetic.f90
)

add_executable(exe src/main.f90)

# Manually link PETSc
target_include_directories(exe PRIVATE ${PETSC_INCLUDES})
target_link_libraries(exe PRIVATE src ${PETSC_LIBRARIES} LAPACK::LAPACK BLAS::BLAS)
#target_include_directories(exe PRIVATE src)
#target_link_libraries(exe PRIVATE src LAPACK::LAPACK BLAS::BLAS)