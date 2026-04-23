#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

module m_petsc
    use m_constants
    use m_basis
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
    Mat, allocatable, intent(inout) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
    Vec, allocatable, intent(inout) :: PROD_VEC(:), FixedSrc(:)
    type(t_mesh), intent(in)      :: mesh
    type(t_finite), intent(in)    :: FE
    type(t_quadrature), intent(in):: Quad
    type(t_material), intent(in)  :: mats(:)
    integer, intent(in)           :: n_groups
    logical, intent(in)           :: is_adjoint

    integer :: ee, g_to, g_from, i, j, k, q, mat_id, row, count, n_valid, i_b, n_basis_patch
    real(dp) :: elem_coords(FE%n_basis, 3), dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
    real(dp) :: detJ, dV, u1, u2, v1, v2, w1, w2
    integer  :: row_b, col_b, row_it, col_it, buf_ptr
    real(dp) :: active_Prod(FE%n_basis), active_Src(FE%n_basis)
    PetscInt :: active_idx(FE%n_basis), idx(FE%n_basis)
    
    real(dp), allocatable :: loc_A(:,:,:)
    real(dp), allocatable :: loc_F(:,:,:,:)
    real(dp), allocatable :: loc_S(:,:,:,:)
    real(dp), allocatable :: loc_Prod(:,:)
    real(dp), allocatable :: loc_Src(:,:)
    real(dp) :: FE_N(FE%n_basis), FE_N_mat(FE%n_basis, FE%n_basis)
    PetscScalar :: val_buf_local(FE%n_basis * FE%n_basis)
    real(dp) :: sigma_s_val, nusigf_val, chi_val

    PetscErrorCode :: ierr
    PetscInt, allocatable :: nnz(:)

    allocate(nnz(mesh%n_nodes)); nnz = 0
    do ee = 1, mesh%n_elems
        n_basis_patch = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee) * mesh%n_cp_zeta(ee)
        do i = 1, n_basis_patch
            row = mesh%elems(ee, i)
            if (row > 0) nnz(row) = nnz(row) + n_basis_patch
        end do
    end do
    ! Cap nnz at total nodes and ensure at least some allocation
    do i = 1, mesh%n_nodes; nnz(i) = max(1, min(nnz(i), mesh%n_nodes)); end do

    if (.not. allocated(A_MAT)) allocate(A_MAT(n_groups))
    if (.not. allocated(MAT_F)) allocate(MAT_F(n_groups, n_groups))
    if (.not. allocated(MAT_S)) allocate(MAT_S(n_groups, n_groups))
    
    if (.not. allocated(PROD_VEC)) allocate(PROD_VEC(n_groups))
    if (.not. allocated(FixedSrc)) allocate(FixedSrc(n_groups))

    do g_to = 1, n_groups
        ! Initialize Vectors
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, PROD_VEC(g_to), ierr)
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, FixedSrc(g_to), ierr)

        ! Create Matrices with precise preallocation
        call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, nnz, A_MAT(g_to), ierr)
        call MatSetOption(A_MAT(g_to), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
        call MatSetOption(A_MAT(g_to), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

        do g_from = 1, n_groups
            call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, nnz, MAT_F(g_to, g_from), ierr)
            call MatSetOption(MAT_F(g_to, g_from), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
            call MatSetOption(MAT_F(g_to, g_from), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

            call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, nnz, MAT_S(g_to, g_from), ierr)
            call MatSetOption(MAT_S(g_to, g_from), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
            call MatSetOption(MAT_S(g_to, g_from), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
        end do
    end do

    count = 0
    write(*,'(A)') " [ MATRIX ] :: Starting element-wise assembly..."
    
    ! Pre-allocate local buffers outside the loop for speed and stability
    allocate(loc_A(FE%n_basis, FE%n_basis, n_groups))
    allocate(loc_F(FE%n_basis, FE%n_basis, n_groups, n_groups))
    allocate(loc_S(FE%n_basis, FE%n_basis, n_groups, n_groups))
    allocate(loc_Prod(FE%n_basis, n_groups), loc_Src(FE%n_basis, n_groups))

    do ee = 1, mesh%n_elems
        loc_A = 0.0_dp; loc_F = 0.0_dp; loc_S = 0.0_dp
        loc_Prod = 0.0_dp; loc_Src = 0.0_dp

        mat_id = mesh%mats(ee)
        n_basis_patch = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee) * mesh%n_cp_zeta(ee)
        
        if (n_basis_patch > FE%n_basis) stop "ERROR: n_basis_patch > FE%n_basis. Check mesh definitions."

        idx = -1
        idx(1:n_basis_patch) = mesh%elems(ee, 1:n_basis_patch) - 1

        do i = 1, n_basis_patch
            row = mesh%elems(ee, i)
            if (row > 0) then
                elem_coords(i, :) = mesh%nodes(row, :)
            else
                elem_coords(i, :) = 0.0_dp
            end if
        end do

        do i = 1, mesh%n_knots_xi_patch(ee) - 1
            u1 = mesh%knot_vectors_xi(ee, i); u2 = mesh%knot_vectors_xi(ee, i+1)
            if (abs(u2 - u1) < 1e-10_dp) cycle
            
            do j = 1, mesh%n_knots_eta_patch(ee) - 1
                v1 = mesh%knot_vectors_eta(ee, j); v2 = mesh%knot_vectors_eta(ee, j+1)
                if (abs(v2 - v1) < 1e-10_dp) cycle

                do k = 1, mesh%n_knots_zeta_patch(ee) - 1
                    w1 = mesh%knot_vectors_zeta(ee, k); w2 = mesh%knot_vectors_zeta(ee, k+1)
                    if (abs(w2 - w1) < 1e-10_dp) cycle

                    do q = 1, Quad%NoPoints
                        call GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, &
                                         elem_coords, dN_dx, dN_dy, dN_dz, detJ, FE_N, FE_N_mat)
                        dV = detJ * Quad%W(q)
                        
                        do g_to = 1, n_groups
                            do row_b = 1, n_basis_patch
                                loc_Prod(row_b, g_to) = loc_Prod(row_b, g_to) + merge(mats(mat_id)%NuSigF(g_to),mats(mat_id)%Chi(g_to), .not. is_adjoint) * FE_N(row_b) * dV
                                loc_Src(row_b, g_to)  = loc_Src(row_b, g_to)  + mats(mat_id)%Src(g_to) * FE_N(row_b) * dV
                                
                                do col_b = 1, n_basis_patch
                                    loc_A(row_b, col_b, g_to) = loc_A(row_b, col_b, g_to) + (mats(mat_id)%D(g_to) * (dN_dx(row_b)*dN_dx(col_b) + dN_dy(row_b)*dN_dy(col_b) + dN_dz(row_b)*dN_dz(col_b)) + mats(mat_id)%SigmaR(g_to) * FE_N_mat(row_b,col_b)) * dV

                                    do g_from = 1, n_groups
                                        nusigf_val = merge(mats(mat_id)%NuSigF(g_from), mats(mat_id)%Chi(g_from), .not. is_adjoint)
                                        chi_val    = merge(mats(mat_id)%Chi(g_to), mats(mat_id)%NuSigF(g_to), .not. is_adjoint)
                                        
                                        loc_F(row_b, col_b, g_to, g_from) = loc_F(row_b, col_b, g_to, g_from) + chi_val * nusigf_val * FE_N_mat(row_b,col_b) * dV

                                        if (g_from /= g_to) then
                                            sigma_s_val = merge(mats(mat_id)%SigmaS(g_from, g_to), mats(mat_id)%SigmaS(g_to, g_from), .not. is_adjoint)
                                            loc_S(row_b, col_b, g_to, g_from) = loc_S(row_b, col_b, g_to, g_from) + sigma_s_val * FE_N_mat(row_b,col_b) * dV
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

        ! Filter out padded indices (-1) before sending to PETSc
        n_valid = 0
        do i = 1, n_basis_patch
            if (idx(i) >= 0) then
                n_valid = n_valid + 1
                active_idx(n_valid) = idx(i)
            end if
        end do

        do g_to = 1, n_groups
            i_b = 0
            do i = 1, FE%n_basis
                if (idx(i) < 0) cycle
                i_b = i_b + 1
                active_Prod(i_b) = loc_Prod(merge(i, 1, i <= n_basis_patch), g_to)
                active_Src(i_b)  = loc_Src(merge(i, 1, i <= n_basis_patch), g_to)
            end do

            call VecSetValues(PROD_VEC(g_to), n_valid, active_idx(1:n_valid), active_Prod(1:n_valid), ADD_VALUES, ierr)
            call VecSetValues(FixedSrc(g_to), n_valid, active_idx(1:n_valid), active_Src(1:n_valid), ADD_VALUES, ierr)
            
            buf_ptr = 1
            do col_it = 1, n_basis_patch
                if (idx(col_it) < 0) cycle
                do row_it = 1, n_basis_patch
                    if (idx(row_it) < 0) cycle
                    val_buf_local(buf_ptr) = loc_A(row_it, col_it, g_to)
                    buf_ptr = buf_ptr + 1
                end do
            end do

            call MatSetValues(A_MAT(g_to), n_valid, active_idx(1:n_valid), n_valid, active_idx(1:n_valid), &
                              val_buf_local, ADD_VALUES, ierr)

            do g_from = 1, n_groups
                ! Pack loc_F into val_buf
                buf_ptr = 1
                do col_it = 1, n_basis_patch
                    if (idx(col_it) < 0) cycle
                    do row_it = 1, n_basis_patch
                        if (idx(row_it) < 0) cycle
                        val_buf_local(buf_ptr) = loc_F(row_it, col_it, g_to, g_from)
                        buf_ptr = buf_ptr + 1
                    end do
                end do
                call MatSetValues(MAT_F(g_to, g_from), n_valid, active_idx(1:n_valid), n_valid, active_idx(1:n_valid), &
                              val_buf_local, ADD_VALUES, ierr)
                if (g_from /= g_to) then
                    ! Pack loc_S into val_buf
                    buf_ptr = 1
                    do col_it = 1, n_basis_patch
                        if (idx(col_it) < 0) cycle
                        do row_it = 1, n_basis_patch
                            if (idx(row_it) < 0) cycle
                            val_buf_local(buf_ptr) = loc_S(row_it, col_it, g_to, g_from)
                            buf_ptr = buf_ptr + 1
                        end do
                    end do
                    call MatSetValues(MAT_S(g_to, g_from), n_valid, active_idx(1:n_valid), n_valid, active_idx(1:n_valid), &
                                      val_buf_local, ADD_VALUES, ierr)
                end if
            end do
        end do

        count = count + 1
        if (mod(count, 7500) == 0) write(*,'(A,I0,A,I0,A,F6.2,A)') ">>> Assembled Element ", count, " of ", mesh%n_elems, " (", real(count,dp)/real(mesh%n_elems,dp)*100.0_dp, "%)"
    end do

    deallocate(loc_A, loc_F, loc_S, loc_Prod, loc_Src)

    ! Finalize vectors, but leave matrices unassembled for BCs
    do i = 1, n_groups
        call VecAssemblyBegin(PROD_VEC(i), ierr); call VecAssemblyEnd(PROD_VEC(i), ierr)
        call VecAssemblyBegin(FixedSrc(i), ierr); call VecAssemblyEnd(FixedSrc(i), ierr)
    end do
    deallocate(nnz)
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
