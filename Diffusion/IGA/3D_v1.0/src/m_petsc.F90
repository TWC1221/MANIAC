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

    integer :: ee, g_to, g_from, i, j, k, q, mat_id, ks, total_spans, node_idx, n_basis_patch
    real(dp) :: detJ, dV, u1, u2, v1, v2, w1, w2, sigma_s_val, nusigf_val, chi_val
    real(dp) :: elem_coords(FE%n_basis, 3), dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
    real(dp) :: FE_N(FE%n_basis), FE_N_mat(FE%n_basis, FE%n_basis)
    integer, allocatable :: span_map(:,:) 
    PetscErrorCode :: ierr
    PetscInt, allocatable :: nnz_array(:)
    integer,  allocatable :: node_to_elem_ptr(:), node_to_elem_list(:)

    allocate(nnz_array(mesh%n_nodes), node_to_elem_ptr(mesh%n_nodes + 1))
    nnz_array = 0; node_to_elem_ptr = 0;
    do ee = 1, mesh%n_elems
        n_basis_patch = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee) * mesh%n_cp_zeta(ee)
        do i = 1, n_basis_patch
            node_idx = mesh%elems(ee, i)
            if (node_idx > 0) node_to_elem_ptr(node_idx) = node_to_elem_ptr(node_idx) + 1
        end do
    end do
    j = 1
    do i = 1, mesh%n_nodes
        k = node_to_elem_ptr(i); node_to_elem_ptr(i) = j; j = j + k
    end do
    node_to_elem_ptr(mesh%n_nodes + 1) = j
    allocate(node_to_elem_list(j - 1))
    do ee = 1, mesh%n_elems
        do i = 1, FE%n_basis
            node_idx = mesh%elems(ee, i)
            if (node_idx > 0) then
                node_to_elem_list(node_to_elem_ptr(node_idx) + nnz_array(node_idx)) = ee
                nnz_array(node_idx) = nnz_array(node_idx) + 1
            end if
        end do
    end do
    do i = 1, mesh%n_nodes
        nnz_array(i) = 0
        do j = node_to_elem_ptr(i), node_to_elem_ptr(i+1) - 1
            nnz_array(i) = nnz_array(i) + FE%n_basis
        end do
        nnz_array(i) = min(nnz_array(i), mesh%n_nodes)
    end do
    deallocate(node_to_elem_ptr, node_to_elem_list)

    total_spans = 0
    do ee = 1, mesh%n_elems
        total_spans = total_spans + (mesh%n_knots_xi_patch(ee)-1)*(mesh%n_knots_eta_patch(ee)-1)*(mesh%n_knots_zeta_patch(ee)-1)
    end do

    allocate(span_map(total_spans, 4)) 
    ks = 0
    do ee = 1, mesh%n_elems
        do k = 1, mesh%n_knots_zeta_patch(ee) - 1
            do j = 1, mesh%n_knots_eta_patch(ee) - 1
                do i = 1, mesh%n_knots_xi_patch(ee) - 1
                    ks = ks + 1
                    span_map(ks, :) = [ee, i, j, k]
                end do
            end do
        end do
    end do

    if (.not. allocated(A_MAT)) allocate(A_MAT(n_groups))
    if (.not. allocated(MAT_F)) allocate(MAT_F(n_groups, n_groups))
    if (.not. allocated(MAT_S)) allocate(MAT_S(n_groups, n_groups))
    if (.not. allocated(PROD_VEC)) allocate(PROD_VEC(n_groups))
    if (.not. allocated(FixedSrc)) allocate(FixedSrc(n_groups))

    do g_to = 1, n_groups
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, PROD_VEC(g_to), ierr)
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, FixedSrc(g_to), ierr)
        call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, nnz_array, A_MAT(g_to), ierr)
        call MatSetOption(A_MAT(g_to), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
        do g_from = 1, n_groups
            call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, nnz_array, MAT_F(g_to, g_from), ierr)
            call MatCreateSeqAIJ(PETSC_COMM_SELF, mesh%n_nodes, mesh%n_nodes, PETSC_NULL_INTEGER, nnz_array, MAT_S(g_to, g_from), ierr)
            call MatSetOption(MAT_F(g_to, g_from), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
            call MatSetOption(MAT_S(g_to, g_from), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
        end do
    end do

    write(*,*) " [ MATRIX ] :: Starting assembly loop over ",(total_spans)," spans..."
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(ks, ee, i, j, k, mat_id, u1, u2, v1, v2, w1, w2, q, detJ, dV, ierr, &
    !$OMP         elem_coords, dN_dx, dN_dy, dN_dz, FE_N, FE_N_mat, &
    !$OMP         sigma_s_val, nusigf_val, chi_val, g_to, g_from)
    do ks = 1, total_spans
        block
            PetscInt :: pidx(FE%n_basis)
            real(dp) :: span_A(FE%n_basis, FE%n_basis), span_F(FE%n_basis, FE%n_basis), span_S(FE%n_basis, FE%n_basis)
            real(dp) :: span_Prod(FE%n_basis), span_Src(FE%n_basis)
            
            ee = span_map(ks, 1); i = span_map(ks, 2); j = span_map(ks, 3); k = span_map(ks, 4)
            mat_id = mesh%mats(ee); pidx = mesh%elems(ee, :) - 1
            
            do q = 1, FE%n_basis
                elem_coords(q, :) = merge(mesh%nodes(mesh%elems(ee, q), :), 0.0_dp, mesh%elems(ee, q) > 0)
            end do

            u1 = mesh%knot_vectors_xi(ee, i);   u2 = mesh%knot_vectors_xi(ee, i+1)
            v1 = mesh%knot_vectors_eta(ee, j);  v2 = mesh%knot_vectors_eta(ee, j+1)
            w1 = mesh%knot_vectors_zeta(ee, k); w2 = mesh%knot_vectors_zeta(ee, k+1)
            if (abs(u2-u1)<1e-10_dp .or. abs(v2-v1)<1e-10_dp .or. abs(w2-w1)<1e-10_dp) cycle

            do g_to = 1, n_groups
                span_A = 0.0_dp; span_Prod = 0.0_dp; span_Src = 0.0_dp
                do q = 1, Quad%NoPoints
                    call GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, elem_coords, dN_dx, dN_dy, dN_dz, detJ, FE_N, FE_N_mat)
                    dV = detJ * Quad%W(q)
                    span_Prod = span_Prod + merge(mats(mat_id)%NuSigF(g_to), mats(mat_id)%Chi(g_to), .not. is_adjoint) * FE_N * dV
                    span_Src  = span_Src  + mats(mat_id)%Src(g_to) * FE_N * dV
                    span_A    = span_A    + (mats(mat_id)%D(g_to) * (spread(dN_dx,2,FE%n_basis)*spread(dN_dx,1,FE%n_basis) + &
                                spread(dN_dy,2,FE%n_basis)*spread(dN_dy,1,FE%n_basis) + spread(dN_dz,2,FE%n_basis)*spread(dN_dz,1,FE%n_basis)) + &
                                mats(mat_id)%SigmaR(g_to) * FE_N_mat) * dV
                end do
                
                !$OMP CRITICAL
                call VecSetValues(PROD_VEC(g_to), FE%n_basis, pidx, span_Prod, ADD_VALUES, ierr)
                call VecSetValues(FixedSrc(g_to), FE%n_basis, pidx, span_Src,  ADD_VALUES, ierr)
                call MatSetValues(A_MAT(g_to), FE%n_basis, pidx, FE%n_basis, pidx, reshape(transpose(span_A), [FE%n_basis ** 2]), ADD_VALUES, ierr)
                !$OMP END CRITICAL

                do g_from = 1, n_groups
                    span_F = 0.0_dp; span_S = 0.0_dp
                    nusigf_val = merge(mats(mat_id)%NuSigF(g_from), mats(mat_id)%Chi(g_from), .not. is_adjoint)
                    chi_val    = merge(mats(mat_id)%Chi(g_to), mats(mat_id)%NuSigF(g_to), .not. is_adjoint)
                    sigma_s_val = merge(mats(mat_id)%SigmaS(g_from, g_to), mats(mat_id)%SigmaS(g_to, g_from), .not. is_adjoint)
                    
                    do q = 1, Quad%NoPoints
                        call GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, elem_coords, dN_dx, dN_dy, dN_dz, detJ, FE_N, FE_N_mat)
                        dV = detJ * Quad%W(q)
                        span_F = span_F + chi_val * nusigf_val * FE_N_mat * dV
                        if (g_from /= g_to) span_S = span_S + sigma_s_val * FE_N_mat * dV
                    end do
                    
                    !$OMP CRITICAL
                    call MatSetValues(MAT_F(g_to, g_from), FE%n_basis, pidx, FE%n_basis, pidx, reshape(transpose(span_F), [FE%n_basis ** 2]), ADD_VALUES, ierr)
                    if (g_from /= g_to) call MatSetValues(MAT_S(g_to, g_from), FE%n_basis, pidx, FE%n_basis, pidx, reshape(transpose(span_S), [FE%n_basis ** 2]), ADD_VALUES, ierr)
                    !$OMP END CRITICAL
                end do
            end do
        end block
    end do
    !$OMP END PARALLEL DO

    ! Finalize vectors, but leave matrices unassembled for BCs
    do i = 1, n_groups
        PetscCall(VecAssemblyBegin(PROD_VEC(i), ierr)); PetscCall(VecAssemblyEnd(PROD_VEC(i), ierr))
        PetscCall(VecAssemblyBegin(FixedSrc(i), ierr)); PetscCall(VecAssemblyEnd(FixedSrc(i), ierr))
    end do
    deallocate(nnz_array, span_map)
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
