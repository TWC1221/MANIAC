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

subroutine SetupKSP_PETSc(ksp_solver, A_mat, ksp_choice, pc_choice)
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
end subroutine SetupKSP_PETSc

subroutine ASSEMBLEMultigroupMAT_PETSc(A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, mesh, FE, Quad, mats, n_groups, is_adjoint)
    Mat, allocatable, intent(out) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
    Vec, allocatable, intent(out) :: PROD_VEC(:), FixedSrc(:)
    type(t_mesh), intent(in)      :: mesh
    type(t_finite), intent(in)    :: FE
    type(t_quadrature), intent(in):: Quad
    type(t_material), intent(in)  :: mats(:)
    integer, intent(in)           :: n_groups
    logical, intent(in)           :: is_adjoint

    integer :: ee, g_to, g_from, i, j, q, mat_id, row, col
    integer, allocatable :: nnz(:)
    real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
    real(dp) :: val_a, val_f, val_s, val_p, val_fixed
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

    allocate(A_MAT(n_groups), MAT_F(n_groups, n_groups), MAT_S(n_groups, n_groups), PROD_VEC(n_groups), FixedSrc(n_groups))

    do g_to = 1, n_groups

        call VecCreate(PETSC_COMM_WORLD, FixedSrc(g_to), ierr)
        call VecSetSizes(FixedSrc(g_to), PETSC_DECIDE, mesh%n_nodes, ierr)
        call VecSetFromOptions(FixedSrc(g_to), ierr)
        call VecSet(FixedSrc(g_to), 0.0_dp, ierr)

        call VecCreate(PETSC_COMM_WORLD, PROD_VEC(g_to), ierr)
        call VecSetSizes(PROD_VEC(g_to), PETSC_DECIDE, mesh%n_nodes, ierr)
        call VecSetFromOptions(PROD_VEC(g_to), ierr)
        call VecSet(PROD_VEC(g_to), 0.0_dp, ierr)

        call MatCreate(PETSC_COMM_WORLD, A_MAT(g_to), ierr)
        call MatSetSizes(A_MAT(g_to), PETSC_DECIDE, PETSC_DECIDE, mesh%n_nodes, mesh%n_nodes, ierr)
        call MatSetFromOptions(A_MAT(g_to), ierr)
        call MatSeqAIJSetPreallocation(A_MAT(g_to), PETSC_NULL_INTEGER, nnz, ierr)

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
                    
                    val_p = merge(mats(mat_id)%NuSigF(g_to),mats(mat_id)%Chi(g_to), .not. is_adjoint) * FE%N(q,i) * dV
                    call VecSetValue(PROD_VEC(g_to), row, val_p, ADD_VALUES, ierr)

                    val_fixed = mats(mat_id)%Src(g_to) * FE%N(q,i) * dV
                    if (abs(val_fixed) > 1e-20_dp) call VecSetValue(FixedSrc(g_to), row, val_fixed, ADD_VALUES, ierr)

                    do j = 1, FE%n_basis
                        col = mesh%elems(ee, j) - 1
                        
                        val_a = (mats(mat_id)%D(g_to) * (dN_dx(i)*dN_dx(j) + dN_dy(i)*dN_dy(j)) + mats(mat_id)%SigmaR(g_to) * FE%N_mat(q,i,j)) * dV
                        call MatSetValue(A_MAT(g_to), row, col, val_a, ADD_VALUES, ierr)

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
        call VecAssemblyBegin(PROD_VEC(g_to), ierr)
        call VecAssemblyEnd(PROD_VEC(g_to), ierr)
        call MatAssemblyBegin(A_MAT(g_to), MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(A_MAT(g_to), MAT_FINAL_ASSEMBLY, ierr)
        do g_from = 1, n_groups
            call MatAssemblyBegin(MAT_F(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(MAT_F(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyBegin(MAT_S(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(MAT_S(g_to, g_from), MAT_FINAL_ASSEMBLY, ierr)
        end do
    end do
end subroutine ASSEMBLEMultigroupMAT_PETSc

subroutine EVALMultigroup_PETSc(B, MAT_F, MAT_S, FixedSrc, X_VEC, X_OLD, &
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
end subroutine EVALMultigroup_PETSc

subroutine EVALProd_PETSc(total_prod, X_VEC, PROD_VEC, n_groups)
    real(dp), intent(out) :: total_prod
    Vec, intent(in)      :: X_VEC(:)   
    Vec, intent(in)      :: PROD_VEC(:)
    integer, intent(in)  :: n_groups
    
    integer :: g, ierr
    PetscScalar :: group_dot
    
    total_prod = 0.0_dp
    do g = 1, n_groups
        call VecDot(PROD_VEC(g), X_VEC(g), group_dot, ierr)
        total_prod = total_prod + real(group_dot, dp)
    end do
end subroutine EVALProd_PETSc

end module m_petsc