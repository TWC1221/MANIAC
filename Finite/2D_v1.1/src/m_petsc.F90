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

subroutine setup_petsc_ksp(ksp_solver, A_mat, ksp_choice, pc_choice)
    KSP, intent(inout)      :: ksp_solver
    Mat, intent(inout)      :: A_mat
    integer, intent(in)     :: ksp_choice, pc_choice
    
    PetscErrorCode :: ierr
    PC             :: pc_ctx

    PetscCall(MatAssemblyBegin(A_mat, MAT_FINAL_ASSEMBLY, ierr))
    PetscCall(MatAssemblyEnd(A_mat, MAT_FINAL_ASSEMBLY, ierr))
    PetscCall(MatSetFromOptions(A_mat, ierr))
    PetscCall(KSPSetOperators(ksp_solver, A_mat, A_mat, ierr))
    
    select case(ksp_choice)
        case(SOLVER_KSP_CG)
            PetscCall(KSPSetType(ksp_solver, KSPCG, ierr))
        case(SOLVER_KSP_GMRES)
            PetscCall(KSPSetType(ksp_solver, KSPGMRES, ierr))
        case(SOLVER_KSP_BCGS)
            PetscCall(KSPSetType(ksp_solver, KSPBCGS, ierr))
    end select

    PetscCall(KSPGetPC(ksp_solver, pc_ctx, ierr))
    select case(pc_choice)
        case(PRECON_NONE)
            PetscCall(PCSetType(pc_ctx, PCNONE, ierr))
        case(PRECON_JACOBI)
            PetscCall(PCSetType(pc_ctx, PCJACOBI, ierr))
        case(PRECON_ILU)
            PetscCall(PCSetType(pc_ctx, PCILU, ierr))
        case(PRECON_CHOLESKY)
            PetscCall(PCSetType(pc_ctx, PCICC, ierr))
        case(PRECON_GAMG)
            PetscCall(PCSetType(pc_ctx, PCGAMG, ierr))
    end select

    PetscCall(KSPSetFromOptions(ksp_solver, ierr))
end subroutine setup_petsc_ksp

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

subroutine assemble_source_matrices_petsc(MAT_F, MAT_S, FixedSrc, mesh, FE, Quad, mats, n_groups, is_adjoint)
    Mat, allocatable, intent(out) :: MAT_F(:,:), MAT_S(:,:)
    Vec, allocatable, intent(out) :: FixedSrc(:)
    type(t_mesh), intent(in)      :: mesh
    type(t_finite), intent(in)    :: FE
    type(t_quadrature), intent(in):: Quad
    type(t_material), intent(in)  :: mats(:)
    integer, intent(in)           :: n_groups
    logical, intent(in)           :: is_adjoint

    integer :: ee, g_to, g_from, i, j, q, mat_id, row, col
    integer, allocatable :: nnz(:)
    real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
    real(dp) :: val_f, val_s, val_fixed
    PetscErrorCode :: ierr
    real(dp) :: sigma_s_val, nusigf_val, chi_val

    allocate(nnz(mesh%n_nodes)); nnz = 0
    do ee = 1, mesh%n_elems
        do i = 1, FE%n_basis
            row = mesh%elems(ee, i)
            nnz(row) = nnz(row) + FE%n_basis
        end do
    end do
    do i = 1, mesh%n_nodes; nnz(i) = min(nnz(i), mesh%n_nodes); end do

    allocate(MAT_F(n_groups, n_groups), MAT_S(n_groups, n_groups), FixedSrc(n_groups))

    do g_to = 1, n_groups
        call VecCreate(PETSC_COMM_WORLD, FixedSrc(g_to), ierr)
        call VecSetSizes(FixedSrc(g_to), PETSC_DECIDE, mesh%n_nodes, ierr)
        call VecSetFromOptions(FixedSrc(g_to), ierr)
        call VecSet(FixedSrc(g_to), 0.0_dp, ierr)

        do g_from = 1, n_groups
            call MatCreate(PETSC_COMM_WORLD, MAT_F(g_to, g_from), ierr)
            call MatSetSizes(MAT_F(g_to, g_from), PETSC_DECIDE, PETSC_DECIDE, mesh%n_nodes, mesh%n_nodes, ierr)
            call MatSetFromOptions(MAT_F(g_to, g_from), ierr)
            call MatSeqAIJSetPreallocation(MAT_F(g_to, g_from), PETSC_NULL_INTEGER, nnz, ierr)

            call MatCreate(PETSC_COMM_WORLD, MAT_S(g_to, g_from), ierr)
            call MatSetSizes(MAT_S(g_to, g_from), PETSC_DECIDE, PETSC_DECIDE, mesh%n_nodes, mesh%n_nodes, ierr)
            call MatSetFromOptions(MAT_S(g_to, g_from), ierr)
            call MatSeqAIJSetPreallocation(MAT_S(g_to, g_from), PETSC_NULL_INTEGER, nnz, ierr)
        end do
    end do

    do ee = 1, mesh%n_elems
        mat_id = mesh%mats(ee)
        do i = 1, FE%n_basis
            elem_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
        end do
        
        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
            dV = detJ * Quad%W(q)
            
            do g_to = 1, n_groups
                do i = 1, FE%n_basis
                    row = mesh%elems(ee, i) - 1
                    
                    val_fixed = mats(mat_id)%Src(g_to) * FE%N(q,i) * dV
                    if (abs(val_fixed) > 1e-20_dp) call VecSetValue(FixedSrc(g_to), row, val_fixed, ADD_VALUES, ierr)

                    do j = 1, FE%n_basis
                        col = mesh%elems(ee, j) - 1
                        
                        do g_from = 1, n_groups
                            if (.not. is_adjoint) then
                                nusigf_val = mats(mat_id)%NuSigF(g_from)
                                chi_val    = mats(mat_id)%Chi(g_to)
                            else
                                nusigf_val = mats(mat_id)%Chi(g_from)
                                chi_val    = mats(mat_id)%NuSigF(g_to)
                            end if
                            val_f = chi_val * nusigf_val * FE%N_mat(q,i,j) * dV
                            if (abs(val_f) > 1e-20_dp) call MatSetValue(MAT_F(g_to, g_from), row, col, val_f, ADD_VALUES, ierr)

                            if (g_from /= g_to) then
                                if (.not. is_adjoint) then
                                    sigma_s_val = mats(mat_id)%SigmaS(g_from, g_to)
                                else
                                    sigma_s_val = mats(mat_id)%SigmaS(g_to, g_from)
                                end if
                                val_s = sigma_s_val * FE%N_mat(q,i,j) * dV
                                if (abs(val_s) > 1e-20_dp) call MatSetValue(MAT_S(g_to, g_from), row, col, val_s, ADD_VALUES, ierr)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do

    deallocate(nnz)
    do g_to = 1, n_groups
        call VecAssemblyBegin(FixedSrc(g_to), ierr)
        call VecAssemblyEnd(FixedSrc(g_to), ierr)
        do g_from = 1, n_groups
            call MatAssemblyBegin(MAT_F(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(MAT_F(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyBegin(MAT_S(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(MAT_S(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
        end do
    end do
end subroutine assemble_source_matrices_petsc

subroutine assemble_multigroup_source_petsc(B, MAT_F, MAT_S, FixedSrc, X_VEC, X_OLD, &
                                            n_groups, k_eff, group_idx, is_eigen)
    Vec, intent(inout) :: B
    Mat, intent(in)    :: MAT_F(:,:), MAT_S(:,:)
    Vec, intent(in)    :: FixedSrc(:)
    Vec, intent(in)    :: X_VEC(:), X_OLD(:)
    integer, intent(in):: n_groups, group_idx
    real(dp), intent(in):: k_eff
    logical, intent(in):: is_eigen
    
    integer :: g, ierr
    Vec :: temp_vec

    call VecDuplicate(B, temp_vec, ierr)
    call VecSet(B, 0.0_dp, ierr)

    if (.not. is_eigen) then
        call VecAXPY(B, 1.0_dp, FixedSrc(group_idx), ierr)
    end if

    do g = 1, n_groups
        if (g /= group_idx) then
            call MatMult(MAT_S(group_idx, g), X_VEC(g), temp_vec, ierr)
            call VecAXPY(B, 1.0_dp, temp_vec, ierr)
        end if

        if (is_eigen) then
            call MatMult(MAT_F(group_idx, g), X_OLD(g), temp_vec, ierr)
            call VecAXPY(B, 1.0_dp / k_eff, temp_vec, ierr)
        end if
    end do

    call VecDestroy(temp_vec, ierr)
end subroutine

subroutine calculate_total_production_petsc(total_prod, X_VEC, mesh, FE, Quad, materials, n_groups)
    real(dp), intent(out)          :: total_prod
    Vec, intent(in)                 :: X_VEC(:)
    type(t_mesh), intent(in)        :: mesh
    type(t_finite), intent(in)      :: FE
    type(t_quadrature), intent(in)  :: Quad
    type(t_material), intent(in)    :: materials(:)
    integer, intent(in)             :: n_groups

    integer :: ee, g, q, i, mat_id, ierr
    real(dp) :: detJ, dV, local_prod
    real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis)
    real(dp) :: local_phi_vals(FE%n_basis)
    integer  :: local_node_indices(FE%n_basis)

    type t_ptr_holder
        PetscScalar, pointer :: p(:)
    end type
    type(t_ptr_holder), allocatable :: x_vec_ptrs(:)

    total_prod = 0.0_dp

    allocate(x_vec_ptrs(n_groups))
    do g = 1, n_groups
        call VecGetArrayRead(X_VEC(g), x_vec_ptrs(g)%p, ierr)
    end do

    do ee = 1, mesh%n_elems
        mat_id = mesh%mats(ee)
        local_node_indices = mesh%elems(ee, :)

        do i = 1, FE%n_basis
            elem_coords(i, :) = mesh%nodes(local_node_indices(i), :)
        end do

        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
            dV = abs(detJ) * Quad%W(q)

            local_prod = 0.0_dp
            do g = 1, n_groups
                local_phi_vals = real(x_vec_ptrs(g)%p(local_node_indices), dp)
                local_prod = local_prod + materials(mat_id)%NuSigF(g) * &
                             dot_product(FE%N(q, :), local_phi_vals)
            end do
            total_prod = total_prod + local_prod * dV
        end do
    end do

    do g = 1, n_groups
        call VecRestoreArrayRead(X_VEC(g), x_vec_ptrs(g)%p, ierr)
    end do
    deallocate(x_vec_ptrs)

end subroutine calculate_total_production_petsc

end module m_petsc