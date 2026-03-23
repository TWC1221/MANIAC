module m_petsc
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
use petscsys
use petscvec
use petscmat
use petscksp
implicit none
private
public :: run_petsc
contains

subroutine run_petsc(m, n, v_csr, j_csr, i_csr, x_vals, b_vals)
  implicit none
  ! Dimensions
  PetscInt, intent(in) :: m, n

  ! CSR arrays (caller-provided)
  PetscInt,    intent(in) :: i_csr(:)   ! size = m+1
  PetscInt,    intent(in) :: j_csr(:)   ! size = nnz
  PetscScalar, intent(in) :: v_csr(:)   ! size = nnz

  ! RHS and solution
  PetscScalar, intent(in)  :: b_vals(:) ! size = m
  PetscScalar, intent(out) :: x_vals(:) ! size = n

  ! Local CSR (0-based, per-row sorted)
  PetscInt,    allocatable :: i0(:), j0(:)
  PetscScalar, allocatable :: v0(:)

  ! PETSc objects
  Mat            :: A
  Vec            :: b, x
  KSP            :: ksp
  PetscErrorCode :: ierr
  PetscBool      :: isInit

  ! Indices for vector operations (0-based)
  PetscInt, allocatable :: idx_b(:), idx_x(:)

  ! Locals for sorting
  PetscInt :: r, p, q, k1, k2, keyj, nnz
  PetscScalar :: keyv

  ! --- Prepare CSR (detect base; convert to 0-based; sort per row) ---
  nnz = size(j_csr)

  allocate(i0(size(i_csr)))
  allocate(j0(nnz))
  allocate(v0(nnz))

  if (i_csr(1) == 0) then
    i0 = i_csr
    j0 = j_csr
  else
    i0 = i_csr - 1
    j0 = j_csr - 1
  end if
  v0 = v_csr

  do r = 1, m
    k1 = i0(r) + 1     
    k2 = i0(r+1)
    if (k2 >= k1) then
      do p = k1 + 1, k2
        keyj = j0(p); keyv = v0(p)
        q = p - 1
        do while (q >= k1.and. j0(q) > keyj)
          j0(q+1) = j0(q)
          v0(q+1) = v0(q)
          q = q - 1
        end do
        j0(q+1) = keyj
        v0(q+1) = keyv
      end do
    end if
  end do

  PetscCall(PetscInitialized(isInit, ierr))
  if (.not. isInit) then
    PetscCall(PetscInitialize(PETSC_NULL_CHARACTER, ierr))
  end if

  ! --- Create Matrix from CSR ---
  PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, m, n, i0, j0, v0, A, ierr))

  ! --- Setup Vectors and Solver ---
  PetscCall(MatCreateVecs(A, x, b, ierr))

  ! Assemble b correctly: insert m entries with 0-based indices
  allocate(idx_b(m))
  idx_b = [(r-1, r=1,m)]
  PetscCall(VecSetValues(b, m, idx_b, b_vals, INSERT_VALUES, ierr))
  PetscCall(VecAssemblyBegin(b, ierr))
  PetscCall(VecAssemblyEnd(b, ierr))

  PetscCall(KSPCreate(PETSC_COMM_SELF, ksp, ierr))
  PetscCall(KSPSetOperators(ksp, A, A, ierr))
  PetscCall(KSPSetFromOptions(ksp, ierr))

  ! --- Solve Ax = b ---
  PetscCall(KSPSolve(ksp, b, x, ierr))

  ! --- Copy solution back to x_vals (no F90 wrappers) ---
  allocate(idx_x(n))
  idx_x = [(r-1, r=1,n)]
  PetscCall(VecGetValues(x, n, idx_x, x_vals, ierr))

  ! Optional view
  ! PetscCall(VecView(x, PETSC_VIEWER_STDOUT_SELF, ierr))

  ! --- Cleanup (do NOT PetscFinalize here) ---
  PetscCall(KSPDestroy(ksp, ierr))
  PetscCall(VecDestroy(x, ierr))
  PetscCall(VecDestroy(b, ierr))
  PetscCall(MatDestroy(A, ierr))

  deallocate(i0, j0, v0, idx_b, idx_x)
end subroutine run_petsc

end module