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

subroutine assemble_multigroup_source_petsc(B, mesh, FE, Quad, mats, X_VEC, X_OLD, &
                                            n_groups, k_eff, group_idx, &
                                            is_eigen, S_ext)
    Vec, intent(inout)                :: B           
    type(t_mesh), intent(in)          :: mesh
    type(t_finite), intent(in)        :: FE
    type(t_quadrature), intent(in)    :: Quad
    type(t_material), intent(in)      :: mats(:)
    Vec, intent(in)                   :: X_VEC(:), X_OLD(:)    
    integer, intent(in)               :: n_groups
    real(dp), intent(in)              :: k_eff
    integer, intent(in)               :: group_idx   
    logical, intent(in)               :: is_eigen
    real(dp), intent(in)              :: S_ext(:,:) 
    
    integer  :: e, q, i, g, ierr, row, mat_id
    real(dp) :: detJ, weight, phi_g_at_q, total_src, total_fission_src, scattering_src
    real(dp) :: local_B(FE%n_basis)
    real(dp) :: elem_coords(FE%n_basis, 2)
    real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis)
    
    call VecSet(B, 0.0_dp, ierr)

    do e = 1, mesh%n_elems
        local_B = 0.0_dp
        mat_id = mesh%mats(e)
        
        do i = 1, FE%n_basis
            elem_coords(i,:) = mesh%nodes(mesh%elems(e,i), :)
        end do

        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
            weight = Quad%W(q) * detJ
            
            total_fission_src = 0.0_dp
            scattering_src    = 0.0_dp
            total_src = 0.0_dp

            do g = 1, n_groups
                call get_local_phi_at_q(X_OLD(g), mesh%elems(e,:), FE%N(q,:), phi_g_at_q)
                total_fission_src = total_fission_src + mats(mat_id)%NuSigF(g) * phi_g_at_q
                if (g /= group_idx) then
                    call get_local_phi_at_q(X_VEC(g), mesh%elems(e,:), FE%N(q,:), phi_g_at_q)
                    scattering_src = scattering_src + mats(mat_id)%SigmaS(group_idx, g) * phi_g_at_q
                end if
            end do

            if (is_eigen) then
                total_src = (mats(mat_id)%Chi(group_idx) / k_eff) * total_fission_src + scattering_src
            else
                total_src = S_ext(e, group_idx) + scattering_src
            end if

            do i = 1, FE%n_basis
                local_B(i) = local_B(i) + FE%N(q,i) * total_src * weight
            end do
        end do

        do i = 1, FE%n_basis
            row = mesh%elems(e,i) - 1
            call VecSetValue(B, row, local_B(i), ADD_VALUES, ierr)
        end do
    end do

    call VecAssemblyBegin(B, ierr)
    call VecAssemblyEnd(B, ierr)
end subroutine

subroutine calculate_total_production_petsc(total_prod, X_VEC, mesh, FE, Quad, materials, n_groups)
    real(dp), intent(out) :: total_prod
    Vec,      intent(in)  :: X_VEC(:)
    type(t_mesh),      intent(in)  :: mesh
    type(t_finite),    intent(in)  :: FE
    type(t_quadrature), intent(in)  :: Quad
    type(t_material),  intent(in)  :: materials(:)
    integer,           intent(in)  :: n_groups

    integer :: ee, i, q, g, mat_id, ierr
    real(dp) :: detJ, dV, phi_at_q, fission_src_at_q
    real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis)
    PetscScalar, pointer :: p_phi(:)
    PetscScalar, pointer :: p_phi_zero(:)

    total_prod = 0.0_dp

    do ee = 1, mesh%n_elems
        mat_id = mesh%mats(ee)
        
        do i = 1, FE%n_basis
            elem_coords(i,:) = mesh%nodes(mesh%elems(ee,i), :)
        end do
        
        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
            dV = detJ * Quad%W(q)
            
            fission_src_at_q = 0.0_dp
            
            do g = 1, n_groups
                call VecGetArrayRead(X_VEC(g), p_phi, ierr)
                p_phi_zero(0:) => p_phi ! Fix 0-based indexing
                
                phi_at_q = 0.0_dp
                do i = 1, FE%n_basis
                    phi_at_q = phi_at_q + p_phi_zero(mesh%elems(ee,i)-1) * FE%N(q,i)
                end do
                
                fission_src_at_q = fission_src_at_q + materials(mat_id)%nuSigF(g) * phi_at_q
                call VecRestoreArrayRead(X_VEC(g), p_phi, ierr)
            end do
            
            total_prod = total_prod + fission_src_at_q * dV
        end do
    end do
end subroutine

subroutine get_local_phi_at_q(X, nodes, basis_q, val)
    Vec, intent(in) :: X
    integer, intent(in) :: nodes(:)
    real(dp), intent(in) :: basis_q(:)
    real(dp), intent(out) :: val
    PetscScalar, pointer :: pX(:)
    integer :: i, ierr

    val = 0.0_dp
    call VecGetArrayRead(X, pX, ierr)
    
    block
        PetscScalar, pointer :: pX_zero(:)
        pX_zero(0:) => pX 
        
        do i = 1, size(nodes)
            val = val + pX_zero(nodes(i)-1) * basis_q(i)
        end do
    end block
    
    call VecRestoreArrayRead(X, pX, ierr)
end subroutine

end module m_petsc