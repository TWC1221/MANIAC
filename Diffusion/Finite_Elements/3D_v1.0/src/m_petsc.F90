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

subroutine ASSEMBLEMultigroupMAT_PETSc(A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, mesh, FE, Quad, mats, n_groups, is_adjoint, printout)
    Mat, allocatable, intent(out) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
    Vec, allocatable, intent(out) :: PROD_VEC(:), FixedSrc(:)
    type(t_mesh), intent(in)      :: mesh
    type(t_finite), intent(in)    :: FE
    type(t_quadrature), intent(in):: Quad
    type(t_material), intent(in)  :: mats(:)
    integer, intent(in)           :: n_groups
    logical, intent(in)           :: is_adjoint
    logical, intent(in), optional :: printout

    PetscInt, allocatable :: nnz_array(:)
    integer, allocatable :: node_to_elem_ptr(:), node_to_elem_list(:)
    integer :: node_idx, elem_idx, local_node, neighbor_node

    integer :: ee, g_to, g_from, i, j, q, mat_id, count
    PetscInt :: ierr
    real(dp) :: elem_coords(FE%n_basis, 3), dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis), detJ, dV
    real(dp) :: sigma_s_val, nusigf_val, chi_val    
    logical :: verbose

    verbose = .false.; if (present(printout)) verbose = printout

    allocate(nnz_array(mesh%n_nodes))
    nnz_array = 0

    allocate(node_to_elem_ptr(mesh%n_nodes + 1))
    node_to_elem_ptr = 0
    do ee = 1, mesh%n_elems
        do i = 1, FE%n_basis
            node_idx = mesh%elems(ee, i)
            node_to_elem_ptr(node_idx) = node_to_elem_ptr(node_idx) + 1
        end do
    end do

    count = 1
    do i = 1, mesh%n_nodes
        j = node_to_elem_ptr(i)
        node_to_elem_ptr(i) = count
        count = count + j
    end do
    node_to_elem_ptr(mesh%n_nodes + 1) = count
    allocate(node_to_elem_list(count - 1))

    do ee = 1, mesh%n_elems
        do i = 1, FE%n_basis
            node_idx = mesh%elems(ee, i)
            node_to_elem_list(node_to_elem_ptr(node_idx) + nnz_array(node_idx)) = ee
            nnz_array(node_idx) = nnz_array(node_idx) + 1
        end do
    end do

    do i = 1, mesh%n_nodes
        nnz_array(i) = 0
        do j = node_to_elem_ptr(i), node_to_elem_ptr(i+1) - 1
            ee = node_to_elem_list(j)
            nnz_array(i) = nnz_array(i) + FE%n_basis
        end do
        nnz_array(i) = min(nnz_array(i), mesh%n_nodes)
    end do

    deallocate(node_to_elem_ptr, node_to_elem_list)

    allocate(A_MAT(n_groups), MAT_F(n_groups, n_groups), MAT_S(n_groups, n_groups), PROD_VEC(n_groups), FixedSrc(n_groups))

    do g_to = 1, n_groups
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, FixedSrc(g_to), ierr)
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, PROD_VEC(g_to), ierr)

        call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, &
                             nnz_array, A_MAT(g_to), ierr)
        call MatSetOption(A_MAT(g_to), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

        do g_from = 1, n_groups
            call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, &
                                 nnz_array, MAT_F(g_to, g_from), ierr)
            call MatSetOption(MAT_F(g_to, g_from), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

            call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, &
                                 nnz_array, MAT_S(g_to, g_from), ierr)
            call MatSetOption(MAT_S(g_to, g_from), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
        end do
    end do

    count = 0
    if (verbose) write(*,'(A)') " [ MATRIX ] :: Starting element-wise assembly..."
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(ee, mat_id, i, j, q, elem_coords, dN_dx, dN_dy, dN_dz, detJ, dV, &
    !$OMP         sigma_s_val, nusigf_val, chi_val, g_to, g_from)
    do ee = 1, mesh%n_elems
        block
        real(dp) :: loc_A(FE%n_basis, FE%n_basis, n_groups)
        real(dp) :: loc_F(FE%n_basis, FE%n_basis, n_groups, n_groups)
        real(dp) :: loc_S(FE%n_basis, FE%n_basis, n_groups, n_groups)
        real(dp) :: loc_Prod(FE%n_basis, n_groups)
        real(dp) :: loc_Src(FE%n_basis, n_groups)
        PetscInt :: idx(FE%n_basis)

        loc_A = 0.0_dp; loc_F = 0.0_dp; loc_S = 0.0_dp
        loc_Prod = 0.0_dp; loc_Src = 0.0_dp

        mat_id = mesh%mats(ee)
        idx = mesh%elems(ee, :) - 1

        do i = 1, FE%n_basis
            elem_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
        end do
        
        do q = 1, Quad%NoPoints
            call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, dN_dz, detJ)

            dV = detJ * Quad%W(q)
            
            do g_to = 1, n_groups
                do i = 1, FE%n_basis
                    loc_Prod(i, g_to) = loc_Prod(i, g_to) + merge(mats(mat_id)%NuSigF(g_to),mats(mat_id)%Chi(g_to), .not. is_adjoint) * FE%N(q,i) * dV
                    loc_Src(i, g_to)  = loc_Src(i, g_to)  + mats(mat_id)%Src(g_to) * FE%N(q,i) * dV
                    
                    do j = 1, FE%n_basis
                        loc_A(i, j, g_to) = loc_A(i, j, g_to) + (mats(mat_id)%D(g_to) * (dN_dx(i)*dN_dx(j) + dN_dy(i)*dN_dy(j) + dN_dz(i)*dN_dz(j)) + mats(mat_id)%SigmaR(g_to) * FE%N_mat(q,i,j)) * dV

                        do g_from = 1, n_groups
                            if (.not. is_adjoint) then
                                nusigf_val = mats(mat_id)%NuSigF(g_from)
                                chi_val    = mats(mat_id)%Chi(g_to)
                            else
                                nusigf_val = mats(mat_id)%Chi(g_from)
                                chi_val    = mats(mat_id)%NuSigF(g_to)
                            end if
                            loc_F(i, j, g_to, g_from) = loc_F(i, j, g_to, g_from) + chi_val * nusigf_val * FE%N_mat(q,i,j) * dV

                            if (g_from /= g_to) then
                                if (.not. is_adjoint) then
                                    sigma_s_val = mats(mat_id)%SigmaS(g_from, g_to)
                                else
                                    sigma_s_val = mats(mat_id)%SigmaS(g_to, g_from)
                                end if
                                loc_S(i, j, g_to, g_from) = loc_S(i, j, g_to, g_from) + sigma_s_val * FE%N_mat(q,i,j) * dV
                            end if
                        end do
                    end do
                end do
            end do
        end do

        !$OMP CRITICAL
        do g_to = 1, n_groups
            call VecSetValues(PROD_VEC(g_to), FE%n_basis, idx, loc_Prod(:, g_to), ADD_VALUES, ierr)
            call VecSetValues(FixedSrc(g_to), FE%n_basis, idx, loc_Src(:, g_to), ADD_VALUES, ierr)
            call MatSetValues(A_MAT(g_to), FE%n_basis, idx, FE%n_basis, idx, reshape(transpose(loc_A(:,:,g_to)), [FE%n_basis**2]), ADD_VALUES, ierr)
            do g_from = 1, n_groups
                call MatSetValues(MAT_F(g_to, g_from), FE%n_basis, idx, FE%n_basis, idx, reshape(transpose(loc_F(:,:,g_to, g_from)), [FE%n_basis**2]), ADD_VALUES, ierr)
                if (g_from /= g_to) then
                    call MatSetValues(MAT_S(g_to, g_from), FE%n_basis, idx, FE%n_basis, idx, reshape(transpose(loc_S(:,:,g_to, g_from)), [FE%n_basis**2]), ADD_VALUES, ierr)
                end if
            end do
        end do
        count = count + 1
        if (verbose .and. mod(count, 7500) == 0) write(*,'(A,I0,A,I0,A,F6.2,A)') ">>> Assembled Element ", count, " of ", mesh%n_elems, " (", real(count,dp)/real(mesh%n_elems,dp)*100.0_dp, "%)"
        !$OMP END CRITICAL
        end block
    end do

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

    if (verbose) write(*,'(A)') " [ MATRIX ] :: Assembly finished."
end subroutine ASSEMBLEMultigroupMAT_PETSc

subroutine EVALMultigroup_PETSc(B, MAT_F, MAT_S, FixedSrc, X_VEC, X_OLD, &
                                            n_groups, k_eff, group_idx, temp_vec, is_eigen)
    Vec, intent(inout) :: B, temp_vec
    Mat, intent(in)    :: MAT_F(:,:), MAT_S(:,:)
    Vec, intent(in)    :: FixedSrc(:)
    Vec, intent(in)    :: X_VEC(:), X_OLD(:)
    integer, intent(in):: n_groups, group_idx
    real(dp), intent(in):: k_eff
    logical, intent(in):: is_eigen
    integer :: g, ierr

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