#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

module m_petsc
    use m_constants
    use m_finite_elements
    use m_material
    use m_quadrature
    use m_types
    
    use petscsys
    use petscvec
    use petscmat
    use petscksp
    
    implicit none

contains

! subroutine create_ksp_solver(ksp, XPETSc, BPETSc, A_MAT, NGNODE, NGRP, PRECOND)
!     KSP, dimension(:), allocatable, intent(inout) :: ksp      
!     Vec, dimension(:), allocatable, intent(inout) :: XPETSc   
!     Vec, dimension(:), allocatable, intent(inout) :: BPETSc   
!     Mat, dimension(:), intent(in)                :: A_MAT    
!     integer, intent(in) :: NGNODE, NGRP, PRECOND  

!     integer         :: GINDX
!     PC              :: pc 
!     PetscErrorCode  :: ierr

!     if (.not. allocated(ksp))    allocate(ksp(NGRP))
!     if (.not. allocated(XPETSc)) allocate(XPETSc(NGRP))
!     if (.not. allocated(BPETSc)) allocate(BPETSc(NGRP))

!     do GINDX = 1, NGRP
!         PetscCall(VecCreate(PETSC_COMM_WORLD, BPETSc(GINDX), ierr))
!         PetscCall(VecSetSizes(BPETSc(GINDX), PETSC_DECIDE, NGNODE, ierr))
!         PetscCall(VecSetFromOptions(BPETSc(GINDX), ierr))

!         PetscCall(VecCreate(PETSC_COMM_WORLD, XPETSc(GINDX), ierr))
!         PetscCall(VecSetSizes(XPETSc(GINDX), PETSC_DECIDE, NGNODE, ierr))
!         PetscCall(VecSetFromOptions(XPETSc(GINDX), ierr))

!         PetscCall(KSPCreate(PETSC_COMM_WORLD, ksp(GINDX), ierr))
!         PetscCall(KSPSetType(ksp(GINDX), KSPCG, ierr))
!         PetscCall(KSPGetPC(ksp(GINDX), pc, ierr))

!         select case(PRECOND)
!             case(0); PetscCall(PCSetType(pc, PCNONE, ierr))
!             case(1); PetscCall(PCSetType(pc, PCJACOBI, ierr))
!             case(2); PetscCall(PCSetType(pc, PCILU, ierr))
!             case(3); PetscCall(PCSetType(pc, PCCHOLESKY, ierr))
!             case(4); PetscCall(PCSetType(pc, PCGAMG, ierr))
!         end select

!         PetscCall(KSPSetOperators(ksp(GINDX), A_MAT(GINDX), A_MAT(GINDX), ierr))
!         PetscCall(KSPSetFromOptions(ksp(GINDX), ierr))
!         PetscCall(KSPSetUp(ksp(GINDX), ierr))
!     end do
! end subroutine create_ksp_solver

subroutine assemble_petsc_matrix(MAT_DATA, A_MAT, mesh, FE, Quad, NGRP, MATIDS)
    type(t_material),   intent(in)  :: MAT_DATA(:)
    type(t_mesh),       intent(in)  :: mesh
    type(t_finite),     intent(in)  :: FE
    type(t_quadrature), intent(in)  :: Quad
    integer,            intent(in)  :: NGRP
    integer,            intent(in)  :: MATIDS(:)
    Mat, allocatable,   intent(out) :: A_MAT(:)
    
    integer :: ee, g, i, j, q, mat_id
    PetscInt :: row, col
    integer, allocatable :: nnz(:)
    real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
    PetscScalar :: val
    PetscErrorCode :: ierr

    allocate(nnz(mesh%n_nodes)); nnz = 0
    do ee = 1, mesh%n_elems
        do i = 1, FE%n_basis
            row = mesh%elems(ee, i)
            nnz(row) = nnz(row) + FE%n_basis
        end do
    end do
    do i = 1, mesh%n_nodes; nnz(i) = min(nnz(i), mesh%n_nodes); end do

    if (allocated(A_MAT)) deallocate(A_MAT)
    allocate(A_MAT(NGRP))

    do g = 1, NGRP
        PetscCall(MatCreate(PETSC_COMM_WORLD, A_MAT(g), ierr))
        PetscCall(MatSetSizes(A_MAT(g), PETSC_DECIDE, PETSC_DECIDE, mesh%n_nodes, mesh%n_nodes, ierr))
        PetscCall(MatSetFromOptions(A_MAT(g), ierr))
        PetscCall(MatSeqAIJSetPreallocation(A_MAT(g), PETSC_NULL_INTEGER, nnz, ierr))
    end do

    do ee = 1, mesh%n_elems
        mat_id = MATIDS(ee)

        do i = 1, FE%n_basis
            elem_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
        end do
        
        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
            dV = detJ * Quad%W(q)
            
            do g = 1, NGRP
               do i = 1, FE%n_basis
                    row = mesh%elems(ee, i) - 1 
                    do j = 1, FE%n_basis
                        col = mesh%elems(ee, j) - 1
                        
                        val = (MAT_DATA(mat_id)%D(g) * (dN_dx(i)*dN_dx(j) + dN_dy(i)*dN_dy(j)) + &
                            MAT_DATA(mat_id)%SigmaR(g) * FE%N_mat(q,i,j)) * dV
                            
                        PetscCall(MatSetValue(A_MAT(g), row, col, val, ADD_VALUES, ierr))
                    end do
                end do
            end do
        end do
    end do

    do g = 1, NGRP
        PetscCall(MatAssemblyBegin(A_MAT(g), MAT_FINAL_ASSEMBLY, ierr))
        PetscCall(MatAssemblyEnd(A_MAT(g), MAT_FINAL_ASSEMBLY, ierr))
    end do
    deallocate(nnz)
end subroutine assemble_petsc_matrix

subroutine assemble_petsc_source_vec(B_VEC, mesh, FE, Quad, NGRP, S_ext)
    type(t_mesh),       intent(in)  :: mesh
    type(t_finite),     intent(in)  :: FE
    type(t_quadrature), intent(in)  :: Quad
    integer,            intent(in)  :: NGRP
    real(dp),           intent(in)  :: S_ext(:,:) 
    Vec, allocatable,   intent(out) :: B_VEC(:)

    integer  :: ee, g, i, q
    PetscInt :: row
    real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
    PetscScalar :: local_val
    PetscErrorCode :: ierr

    if (allocated(B_VEC)) deallocate(B_VEC)
    allocate(B_VEC(NGRP))

    do g = 1, NGRP
        PetscCall(VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, B_VEC(g), ierr))
        PetscCall(VecSetFromOptions(B_VEC(g), ierr))
    end do

    do ee = 1, mesh%n_elems
        do i = 1, FE%n_basis
            elem_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
        end do

        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
            dV = detJ * Quad%W(q)
            
            do g = 1, NGRP
                do i = 1, FE%n_basis
                    row = mesh%elems(ee, i) - 1 
                    local_val = S_ext(ee, g) * FE%N(q, i) * dV
                    PetscCall(VecSetValue(B_VEC(g), row, local_val, ADD_VALUES, ierr))
                end do
            end do
        end do
    end do

    do g = 1, NGRP
        PetscCall(VecAssemblyBegin(B_VEC(g), ierr))
        PetscCall(VecAssemblyEnd(B_VEC(g), ierr))
    end do
end subroutine assemble_petsc_source_vec

subroutine assemble_multigroup_fission_vec(Q_PETSC, mesh, FE, Quad, materials, X_VEC, n_groups, k_eff)
    Vec, intent(inout)              :: Q_PETSC
    type(t_mesh), intent(in)        :: mesh
    type(t_finite), intent(in)      :: FE
    type(t_quadrature), intent(in)  :: Quad
    type(t_material), intent(in)    :: materials(:)
    Vec, intent(in)                 :: X_VEC(:)
    integer, intent(in)             :: n_groups
    real(dp), intent(in)            :: k_eff

    integer :: ee, g, i, q, mat_id, ierr
    real(dp) :: detJ, dV, local_phi_fission, inv_k
    real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis)
    real(dp) :: phi_at_nodes(FE%n_basis)
    integer  :: idx(FE%n_basis)
    PetscScalar, allocatable :: local_Q(:)
    PetscScalar, pointer     :: p_phi(:) ! Pointer for direct memory access

    inv_k = 1.0_dp / k_eff
    call VecZeroEntries(Q_PETSC, ierr)
    allocate(local_Q(FE%n_basis))

    do ee = 1, mesh%n_elems
        mat_id = mesh%mats(ee)
        idx = mesh%elems(ee, :) - 1
        
        do i = 1, FE%n_basis
            elem_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
        end do

        local_Q = 0.0_dp

        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
            dV = detJ * Quad%W(q)
            
            local_phi_fission = 0.0_dp
            do g = 1, n_groups
                call VecGetArrayRead(X_VEC(g), p_phi, ierr)
                do i = 1, FE%n_basis
                    phi_at_nodes(i) = p_phi(idx(i) + 1)
                end do
                call VecRestoreArrayRead(X_VEC(g), p_phi, ierr)
                
                local_phi_fission = local_phi_fission + &
                    materials(mat_id)%nuSigF(g) * dot_product(FE%N(q, :), phi_at_nodes)
            end do

            do i = 1, FE%n_basis
                local_Q(i) = local_Q(i) + inv_k * local_phi_fission * FE%N(q, i) * dV
            end do
        end do

        call VecSetValues(Q_PETSC, FE%n_basis, idx, local_Q, ADD_VALUES, ierr)
    end do

    call VecAssemblyBegin(Q_PETSC, ierr)
    call VecAssemblyEnd(Q_PETSC, ierr)
    deallocate(local_Q)
end subroutine assemble_multigroup_fission_vec

subroutine normalize_flux(X_VEC, n_groups)
    Vec, intent(inout)  :: X_VEC(:)
    integer, intent(in) :: n_groups
    integer             :: g, ierr
    real(dp)            :: group_norm, total_norm_sq, global_norm

    total_norm_sq = 0.0_dp

    do g = 1, n_groups
        call VecNorm(X_VEC(g), NORM_2, group_norm, ierr)
        total_norm_sq = total_norm_sq + group_norm**2
    end do
    
    global_norm = sqrt(total_norm_sq)

    if (global_norm > 1e-18) then
        do g = 1, n_groups
            call VecScale(X_VEC(g), 1.0_dp / global_norm, ierr)
        end do
    end if
end subroutine normalize_flux
end module m_petsc