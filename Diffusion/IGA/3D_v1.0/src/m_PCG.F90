
!------------------------------------------------------------------------!
! Purpose:                                                             -!
!  Carry out sparse matrix multiplication in CSR format                 -!  
!  Carry out Conjugate Gradient method for solving linear systems       -!  
!  Generate Preconditioning matrices M for use in Preconditioned CG     -!  
!  using Incomplete Cholesky factorization (IC(0)) in CSR format        -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 16/02/25     T. Charlton    Implemented 2D Application                -!
!------------------------------------------------------------------------! 

module m_PCG
use m_constants
use m_basis
use m_material
use m_quadrature
use m_types
implicit none

contains 

subroutine assemble_PCG_matrix(MAT_DATA, A_MAT_PCG, mesh, FE, Quad, NGRP, MATIDS)
    type(t_material),       intent(in)  :: MAT_DATA(:)
    type(t_mesh),           intent(in)  :: mesh
    type(t_finite),         intent(in)  :: FE
    type(t_quadrature),     intent(in)  :: Quad
    integer,                intent(in)  :: NGRP
    integer,                intent(in)  :: MATIDS(:)
    type(t_PCG_CSR), allocatable,  intent(out) :: A_MAT_PCG(:)
    
    integer :: ee, g, i, j, q, mat_id, row, row_b, col, col_b, n_basis_patch, node_id
    real(dp) :: FE_N(FE%n_basis), FE_N_mat(FE%n_basis, FE%n_basis)
    integer, allocatable :: nnz(:) ; integer :: k
    real(dp) :: u1, u2, v1, v2, w1, w2
    real(dp) :: elem_coords(FE%n_basis, 3), dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis), detJ, dV, val

    allocate(nnz(mesh%n_nodes)); nnz = 0
    do ee = 1, mesh%n_elems
        do i = 1, FE%n_basis
            row = mesh%elems(ee, i)
            nnz(row) = nnz(row) + FE%n_basis
        end do
    end do
    do i = 1, mesh%n_nodes; nnz(i) = min(nnz(i), mesh%n_nodes); end do

    if (allocated(A_MAT_PCG)) deallocate(A_MAT_PCG)
    allocate(A_MAT_PCG(NGRP))

    do g = 1, NGRP
        call PCG_MAT_INIT(A_MAT_PCG(g), mesh, nnz)
        A_MAT_PCG(g)%row_ptr(1) = 1
        do i = 1, mesh%n_nodes
            A_MAT_PCG(g)%row_ptr(i+1) = A_MAT_PCG(g)%row_ptr(i) + nnz(i)
        end do
    end do

    do ee = 1, mesh%n_elems
        mat_id = MATIDS(ee)
        n_basis_patch = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee) * mesh%n_cp_zeta(ee)

        do i = 1, n_basis_patch
            node_id = mesh%elems(ee, i)
            elem_coords(i, 1:3) = merge(mesh%nodes(node_id, 1:3), 0.0_dp, node_id > 0)
        end do
        
        do i = 1, size(mesh%knot_vectors_xi, 2) - 1
            u1 = mesh%knot_vectors_xi(ee, i); u2 = mesh%knot_vectors_xi(ee, i+1)
            if (abs(u2 - u1) < 1e-10_dp) cycle
            do j = 1, size(mesh%knot_vectors_eta, 2) - 1
                v1 = mesh%knot_vectors_eta(ee, j); v2 = mesh%knot_vectors_eta(ee, j+1)
                if (abs(v2 - v1) < 1e-10_dp) cycle

                do k = 1, size(mesh%knot_vectors_zeta, 2) - 1
                    w1 = mesh%knot_vectors_zeta(ee, k); w2 = mesh%knot_vectors_zeta(ee, k+1)
                    if (abs(w2 - w1) < 1e-10_dp) cycle

                    do q = 1, Quad%NoPoints
                        call GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, &
                                         elem_coords, dN_dx, dN_dy, dN_dz, detJ, FE_N, FE_N_mat)
                        dV = detJ * Quad%W(q)
                        
                        do g = 1, NGRP
                            do row_b = 1, n_basis_patch
                                row = mesh%elems(ee, row_b)
                                do col_b = 1, n_basis_patch
                                    col = mesh%elems(ee, col_b)
                                    
                                    val = (MAT_DATA(mat_id)%D(g) * (dN_dx(row_b)*dN_dx(col_b) + dN_dy(row_b)*dN_dy(col_b) + &
                                        dN_dz(row_b)*dN_dz(col_b)) + MAT_DATA(mat_id)%SigmaR(g) * FE_N_mat(row_b,col_b)) * dV
                                        
                                    call PCG_MAT_ALLOCATION(A_MAT_PCG(g), row, col, val)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    do g = 1, NGRP
        call PCG_COMPRESS_CSR(A_MAT_PCG(g), mesh%n_nodes)
        call Sort_CSR(A_MAT_PCG(g))
    end do

    deallocate(nnz)
end subroutine assemble_PCG_matrix

subroutine assemble_source_matrices_pcg(MAT_F, MAT_S, FixedSrc, mesh, FE, Quad, mats, n_groups, is_adjoint)
    type(t_PCG_CSR), allocatable, intent(out) :: MAT_F(:,:), MAT_S(:,:)
    type(t_vec), allocatable, intent(out) :: FixedSrc(:)
    type(t_mesh), intent(in)      :: mesh
    type(t_finite), intent(in)    :: FE
    type(t_quadrature), intent(in):: Quad
    type(t_material), intent(in)  :: mats(:)
    integer, intent(in)           :: n_groups
    logical, intent(in)           :: is_adjoint

    integer :: ee, g_to, g_from, i, j, k, q, mat_id, row, row_b, col, col_b, n_basis_patch
    real(dp) :: FE_N(FE%n_basis), FE_N_mat(FE%n_basis, FE%n_basis)
    integer, allocatable :: nnz(:)
    real(dp) :: elem_coords(FE%n_basis, 3), dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis), detJ, dV
    real(dp) :: val_f, val_s, val_fixed
    real(dp) :: u1, u2, v1, v2, w1, w2
    real(dp) :: sigma_s_val, nusigf_val, chi_val

    allocate(nnz(mesh%n_nodes)); nnz = 0
    do ee = 1, mesh%n_elems
        n_basis_patch = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee) * mesh%n_cp_zeta(ee)
        do i = 1, n_basis_patch
            row = mesh%elems(ee, i)
            if (row > 0) nnz(row) = nnz(row) + n_basis_patch
        end do
    end do
    do i = 1, mesh%n_nodes; nnz(i) = min(nnz(i), mesh%n_nodes); end do

    allocate(MAT_F(n_groups, n_groups), MAT_S(n_groups, n_groups), FixedSrc(n_groups))

    do g_to = 1, n_groups
        call PCG_VEC_INIT(mesh%n_nodes, FixedSrc(g_to))
        do g_from = 1, n_groups
            call PCG_MAT_INIT(MAT_F(g_to, g_from), mesh, nnz)
            call PCG_MAT_INIT(MAT_S(g_to, g_from), mesh, nnz)
            MAT_F(g_to, g_from)%row_ptr(1) = 1
            MAT_S(g_to, g_from)%row_ptr(1) = 1
            do i = 1, mesh%n_nodes
                MAT_F(g_to, g_from)%row_ptr(i+1) = MAT_F(g_to, g_from)%row_ptr(i) + nnz(i)
                MAT_S(g_to, g_from)%row_ptr(i+1) = MAT_S(g_to, g_from)%row_ptr(i) + nnz(i)
            end do
        end do
    end do

    do ee = 1, mesh%n_elems
        mat_id = mesh%mats(ee)
        n_basis_patch = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee) * mesh%n_cp_zeta(ee)
        do i = 1, FE%n_basis
            if (mesh%elems(ee, i) > 0) then
                elem_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
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
                                row = mesh%elems(ee, row_b)
                                val_fixed = mats(mat_id)%Src(g_to) * FE_N(row_b) * dV
                                if (abs(val_fixed) > 1e-20_dp) call PCG_VEC_ALLOCATION(FixedSrc(g_to), row, val_fixed)

                                do col_b = 1, n_basis_patch
                                    col = mesh%elems(ee, col_b)
                                    do g_from = 1, n_groups
                                        nusigf_val = merge(mats(mat_id)%NuSigF(g_from), mats(mat_id)%Chi(g_from), .not. is_adjoint)
                                        chi_val    = merge(mats(mat_id)%Chi(g_to), mats(mat_id)%NuSigF(g_to), .not. is_adjoint)
                                        
                                        val_f = chi_val * nusigf_val * FE_N_mat(row_b,col_b) * dV
                                        if (abs(val_f) > 1e-20_dp) call PCG_MAT_ALLOCATION(MAT_F(g_to, g_from), row, col, val_f)

                                        if (g_from /= g_to) then
                                            sigma_s_val = merge(mats(mat_id)%SigmaS(g_from, g_to), mats(mat_id)%SigmaS(g_to, g_from), .not. is_adjoint)
                                            val_s = sigma_s_val * FE_N_mat(row_b,col_b) * dV
                                            if (abs(val_s) > 1e-20_dp) call PCG_MAT_ALLOCATION(MAT_S(g_to, g_from), row, col, val_s)
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    do g_to = 1, n_groups
        do g_from = 1, n_groups
            call PCG_COMPRESS_CSR(MAT_F(g_to, g_from), mesh%n_nodes)
            call Sort_CSR(MAT_F(g_to, g_from))
            call PCG_COMPRESS_CSR(MAT_S(g_to, g_from), mesh%n_nodes)
            call Sort_CSR(MAT_S(g_to, g_from))
        end do
    end do
    deallocate(nnz)
end subroutine assemble_source_matrices_pcg

subroutine assemble_multigroup_source_pcg(B, MAT_F, MAT_S, FixedSrc, X_VEC, X_OLD, &
                                            n_groups, k_eff, group_idx, is_eigen)
    type(t_vec), intent(inout)        :: B           
    type(t_PCG_CSR), intent(in)           :: MAT_F(:,:), MAT_S(:,:)
    type(t_vec), intent(in)           :: FixedSrc(:)
    type(t_vec), intent(in)           :: X_VEC(:), X_OLD(:)
    integer, intent(in)               :: n_groups, group_idx
    real(dp), intent(in)              :: k_eff
    logical, intent(in)               :: is_eigen
    
    integer :: g

    B%vec = 0.0_dp
    if (.not. is_eigen) B%vec = FixedSrc(group_idx)%vec
    
    do g = 1, n_groups
        if (g /= group_idx) then
            B%vec = B%vec + CSR_dot_product(MAT_S(group_idx, g)%val, MAT_S(group_idx, g)%col, &
                                            MAT_S(group_idx, g)%row_ptr, X_VEC(g)%vec)
        end if
    end do

    do g = 1, n_groups
        B%vec = B%vec + (1.0_dp / k_eff) * CSR_dot_product(MAT_F(group_idx, g)%val, MAT_F(group_idx, g)%col, &
                                                           MAT_F(group_idx, g)%row_ptr, X_OLD(g)%vec)
    end do
end subroutine assemble_multigroup_source_pcg

subroutine PCG_MAT_INIT(A_MAT_PCG, mesh, nnz_array)
    type(t_mesh), intent(in) :: mesh
    type(t_PCG_CSR),  intent(inout) :: A_MAT_PCG
    integer,      intent(in)    :: nnz_array(:)
    integer :: total_nnz

    total_nnz = sum(nnz_array)

    if (allocated(A_MAT_PCG%val)) deallocate(A_MAT_PCG%val)
    if (allocated(A_MAT_PCG%col)) deallocate(A_MAT_PCG%col)
    if (allocated(A_MAT_PCG%row_ptr)) deallocate(A_MAT_PCG%row_ptr)

    allocate(A_MAT_PCG%val(total_nnz))
    allocate(A_MAT_PCG%col(total_nnz))
    allocate(A_MAT_PCG%row_ptr(mesh%n_nodes + 1))
    
    A_MAT_PCG%val = 0.0_dp
    A_MAT_PCG%col = 0
    A_MAT_PCG%row_ptr = 0
end subroutine PCG_MAT_INIT

subroutine PCG_MAT_ALLOCATION(A_MAT_PCG, row, col, val) ! Note: row is not used for lookup in current CSR builder logic, assumed correct by caller loop
    type(t_PCG_CSR), intent(inout) :: A_MAT_PCG
    integer,     intent(in)    :: row, col
    real(dp),    intent(in)    :: val
    integer                    :: k, row_start, row_end

    row_start = A_MAT_PCG%row_ptr(row)
    row_end   = A_MAT_PCG%row_ptr(row+1) - 1

    do k = row_start, row_end
        if (A_MAT_PCG%col(k) == col) then
            A_MAT_PCG%val(k) = A_MAT_PCG%val(k) + val
            return
        end if
        
        if (A_MAT_PCG%col(k) == 0) then
            A_MAT_PCG%col(k) = col
            A_MAT_PCG%val(k) = val
            return
        end if
    end do
    
end subroutine PCG_MAT_ALLOCATION

subroutine PCG_COMPRESS_CSR(A_MAT_PCG, n_nodes)
    type(t_PCG_CSR), intent(inout) :: A_MAT_PCG
    integer,     intent(in)    :: n_nodes
    integer, allocatable       :: new_row_ptr(:)
    integer                    :: i
    integer,  allocatable :: tmp_c(:)
    real(dp), allocatable :: tmp_v(:)

    if (allocated(new_row_ptr)) deallocate(new_row_ptr)
    allocate(new_row_ptr(n_nodes + 1))
    new_row_ptr(1) = 1

    do i = 1, n_nodes
        new_row_ptr(i+1) = new_row_ptr(i) + count(A_MAT_PCG%col(A_MAT_PCG%row_ptr(i) : A_MAT_PCG%row_ptr(i+1)-1) /= 0)
    end do
    
    tmp_v = pack(A_MAT_PCG%val, A_MAT_PCG%col /= 0)
    tmp_c = pack(A_MAT_PCG%col, A_MAT_PCG%col /= 0)
    
    call move_alloc(tmp_v, A_MAT_PCG%val)
    call move_alloc(tmp_c, A_MAT_PCG%col)
    call move_alloc(new_row_ptr, A_MAT_PCG%row_ptr)
end subroutine PCG_COMPRESS_CSR

subroutine Sort_CSR(A)
    type(t_PCG_CSR), intent(inout) :: A
    integer :: i, j, k, row_start, row_end, n_nodes
    integer :: temp_c
    real(dp) :: temp_v
    logical :: swapped

    n_nodes = size(A%row_ptr) - 1

    do i = 1, n_nodes
        row_start = A%row_ptr(i)
        row_end   = A%row_ptr(i+1) - 1
        
        ! Simple Bubble Sort for the row (bandwidth is usually small)
        do j = row_end, row_start + 1, -1
            swapped = .false.
            do k = row_start, j - 1
                if (A%col(k) > A%col(k+1)) then
                    temp_c = A%col(k)
                    A%col(k) = A%col(k+1)
                    A%col(k+1) = temp_c
                    
                    temp_v = A%val(k)
                    A%val(k) = A%val(k+1)
                    A%val(k+1) = temp_v
                    swapped = .true.
                end if
            end do
            if (.not. swapped) exit
        end do
    end do
end subroutine Sort_CSR

subroutine PCG_VEC_INIT(n_nodes, B_STRUCT)
    integer, intent(in) :: n_nodes
    type(t_VEC), intent(inout) :: B_STRUCT

    if (allocated(B_STRUCT%vec)) deallocate(B_STRUCT%vec)
    allocate(B_STRUCT%vec(n_nodes))
    
    B_STRUCT%vec = 0.0_dp
end subroutine PCG_VEC_INIT

subroutine PCG_VEC_ALLOCATION(B_STRUCT, row, local_val)
    type(t_VEC), intent(inout) :: B_STRUCT
    integer,     intent(in)    :: row
    real(dp),    intent(in)    :: local_val

    B_STRUCT%vec(row) = B_STRUCT%vec(row) + local_val
end subroutine PCG_VEC_ALLOCATION

subroutine calculate_total_production_pcg(total_prod, X_PCG, mesh, FE, Quad, materials, n_groups)
    real(dp), intent(out)           :: total_prod
    type(t_vec), intent(in)         :: X_PCG(:)
    type(t_mesh), intent(in)        :: mesh
    type(t_finite), intent(in)      :: FE
    type(t_quadrature), intent(in)  :: Quad
    type(t_material), intent(in)    :: materials(:)
    integer, intent(in)             :: n_groups

    integer  :: ee, g, i, j, k, q, mat_id, n_basis_patch
    integer  :: node_indices(FE%n_basis)
    real(dp) :: detJ, dV, local_phi, fission_rate_at_q
    real(dp) :: FE_N(FE%n_basis)
    real(dp) :: elem_coords(FE%n_basis, 3)
    real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
    real(dp) :: u1, u2, v1, v2, w1, w2

    total_prod = 0.0_dp

    do ee = 1, mesh%n_elems
        mat_id = mesh%mats(ee)
        n_basis_patch = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee) * mesh%n_cp_zeta(ee)
        node_indices = mesh%elems(ee, :)
        
        do i = 1, FE%n_basis
            if (node_indices(i) > 0) elem_coords(i, :) = mesh%nodes(node_indices(i), :)
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
                                         elem_coords, dN_dx, dN_dy, dN_dz, detJ, FE_N)
                        dV = detJ * Quad%W(q)
                        
                        fission_rate_at_q = 0.0_dp
                        do g = 1, n_groups
                            local_phi = dot_product(X_PCG(g)%vec(node_indices(1:n_basis_patch)), FE_N(1:n_basis_patch))
                            fission_rate_at_q = fission_rate_at_q + materials(mat_id)%NuSigF(g) * local_phi
                        end do
                        total_prod = total_prod + fission_rate_at_q * dV
                    end do
                end do
            end do
        end do
    end do
end subroutine calculate_total_production_pcg

subroutine calculate_mesh_volume(total_vol, mat_vols, mesh, FE, Quad)
    real(dp), intent(out)           :: total_vol
    real(dp), intent(out), allocatable :: mat_vols(:)
    type(t_mesh), intent(in)        :: mesh
    type(t_finite), intent(in)      :: FE
    type(t_quadrature), intent(in)  :: Quad

    integer  :: ee, i, j, k, q, node_id, mat_id, max_mat
    real(dp) :: detJ, dV, u1, u2, v1, v2, w1, w2
    real(dp) :: elem_coords(FE%n_basis, 3)
    real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis), FE_N(FE%n_basis)

    total_vol = 0.0_dp
    max_mat = maxval(mesh%mats)
    allocate(mat_vols(max_mat))
    mat_vols = 0.0_dp

    do ee = 1, mesh%n_elems
        mat_id = mesh%mats(ee)
        elem_coords = 0.0_dp
        do i = 1, FE%n_basis
            node_id = mesh%elems(ee, i)
            if (node_id > 0) elem_coords(i, :) = mesh%nodes(node_id, :)
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
                                         elem_coords, dN_dx, dN_dy, dN_dz, detJ, FE_N)
                        dV = detJ * Quad%W(q)
                        total_vol = total_vol + dV
                        if (mat_id > 0) mat_vols(mat_id) = mat_vols(mat_id) + dV
                    end do
                end do
            end do
        end do
    end do
end subroutine calculate_mesh_volume

    !-------------------------
    ! Principle Conjugation Gradient (PCG) Algorithm driver, solving Ax=b
    !-------------------------
    subroutine CGSolve(A_MAT_PCG, x, b, PCG_mode, max_iter)
        implicit none 
        type(t_PCG_CSR), intent(in) :: A_MAT_PCG
        type(t_vec), intent(inout) :: x, b
        real(8) :: alpha, beta
        integer :: ii
        integer, intent(in)  :: PCG_mode, max_iter

        real(8) :: rho_old, rho_new, denom
        real(8), allocatable :: L_AA(:), U_AA(:), diag(:), r(:), d(:), z(:), q(:)
        integer, allocatable :: L_IA(:), L_JA(:), U_IA(:), U_JA(:)
        logical :: nan_detected

        associate(AA => A_MAT_PCG%val, JA => A_MAT_PCG%col, IA => A_MAT_PCG%row_ptr, &
                  xv => x%vec, bv => b%vec)

        if (.not. allocated(r)) allocate(r(size(bv)))
        if (.not. allocated(d)) allocate(d(size(bv)))
        if (.not. allocated(z)) allocate(z(size(bv)))
        if (.not. allocated(q)) allocate(q(size(bv)))

            r = bv - CSR_dot_product(AA, JA, IA, xv)

        select case(PCG_mode)
            case(PRECON_NONE)
                z = r
            case(PRECON_CHOLESKY)
                call Cholesky_CSR(AA, JA, IA, L_AA, L_JA, L_IA)
                z = PCG_Cholesky_CSR(L_AA, L_JA, L_IA, r)
            case(PRECON_ILU)
                call ILU0_CSR(AA, JA, IA, L_AA, L_JA, L_IA, U_AA, U_JA, U_IA)
                z = PCG_ILU_CSR(L_AA, L_JA, L_IA, U_AA, U_JA, U_IA, r)
            case(PRECON_JACOBI)
                call Jacobi_CSR(AA, JA, IA, diag)
                z = r/diag
        end select

        ! Check for NaN in initial z
        nan_detected = .false.
        if (any(isnan(z))) then
            print *, "WARNING: NaN detected in initial preconditioned residual z"
            z = r  
            nan_detected = .true.
        end if

        d = z
        ! rho = r^T z (used to compute alpha and beta)
        rho_old = dot_product(r, z)

        ! if (isnan(rho_old)) then
        !     print *, "WARNING: rho_old is NaN in PCG initialization"
        !     rho_old = 1.0d0
        ! end if

        do ii = 1, max_iter
            q = CSR_dot_product(AA, JA, IA, d)
            denom = dot_product(d, q)
            
            ! if (abs(denom) < 1.0d-20) then
            !     print *, "WARNING: denom too small in PCG iteration", ii, ":", denom
            !     print *, "  rho_old =", rho_old
            !     print *, "  ||d|| =", sqrt(dot_product(d,d))
            !     print *, "  ||q|| =", sqrt(dot_product(q,q))
            !     exit
            ! end if
            
            alpha = rho_old / denom

            if (isnan(alpha)) then
                print *, "ERROR: alpha is NaN in iteration", ii
                print *, "  rho_old =", rho_old
                print *, "  denom =", denom
                print *, "  ||d|| =", sqrt(dot_product(d,d))
                print *, "  ||q|| =", sqrt(dot_product(q,q))
                print *, "  Residual norm =", sqrt(dot_product(r,r))
                exit
            end if

            if (abs(alpha) > 1.0d10) then
                print *, "WARNING: alpha too large in iteration", ii, ":", alpha
                exit
            end if

            x%vec = x%vec + alpha*d
            r = r - alpha*q

            ! Apply preconditioner to new residual
            select case(PCG_mode)
                case(PRECON_NONE)
                    z = r
                case(PRECON_CHOLESKY)
                    z = PCG_Cholesky_CSR(L_AA, L_JA, L_IA, r)
                case(PRECON_ILU)
                    z = PCG_ILU_CSR(L_AA, L_JA, L_IA, U_AA, U_JA, U_IA, r)
                case(PRECON_JACOBI)
                    z = r / diag
            end select

            if (any(isnan(z))) then
                print *, "WARNING: NaN detected in preconditioned residual at iteration", ii
                print *, "  Setting z = r (unpreconditioned)"
                z = r
            end if

            rho_new = dot_product(r, z)
            
            if (isnan(rho_new)) then
                print *, "WARNING: rho_new is NaN in iteration", ii
                exit
            end if

            if (abs(rho_old) < 1.0d-24) then
                !print *, "WARNING: rho_old too small in iteration", ii, ":", rho_old
                exit
            end if
            
            !print*, ii, r(1), r(50), r(2000), r(4000)
            
            beta = rho_new / rho_old
            d = z + beta*d

            if (rho_new < 1.0d-10 * rho_old) exit ! Converged!

            rho_old = rho_new
        end do

        if (any(isnan(x%vec))) then
            print *, "ERROR: NaN detected in final solution x"
        end if

        end associate
    end subroutine

    !-------------------------
    ! CSR Matvec multiplication function
    !-------------------------
    function CSR_dot_product(AA,JA,IA,x) result(Y)
        implicit none
        real(8), intent(in) :: AA(:), x(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8) :: Y(SIZE(x))
        integer :: ii, k1, k2, n_rows
        n_rows = size(IA) - 1

        do ii = 1, n_rows
            k1 = IA(ii)
            k2 = IA(ii+1)-1
            Y(ii) = dot_product(AA(k1:k2),x(JA(k1:k2)))
        end do
        return
    end function

    !-------------------------
    ! Incomplete Cholesky (lower, CSR) Decomposition of CSR matrix A
    ! AA, JA, IA define original A in CSR. 
    ! L_AA, L_JA, L_IA define the lower triangular L. L transformation into Cholesky preconditioning matrix M done in place
    ! Returns L_AA, L_JA, L_IA. The CSR format of ICD(0) matrix M
    !-------------------------
    subroutine Cholesky_CSR(AA, JA, IA, L_AA, L_JA, L_IA)
        real(8), intent(in) :: AA(:)
        integer, intent(in) :: JA(:), IA(:)

        real(8), allocatable, intent(out) :: L_AA(:)
        integer, allocatable, intent(out) :: L_JA(:), L_IA(:)

        real(8), allocatable :: L_D0(:)
        real(8) :: s
        integer :: ii, jj, kk, k1, k2, nnz, jcol, kcol, kk1, kk2, n
        
        n = size(IA) - 1
        nnz = 0

        do ii = 1, n
            k1 = IA(ii)
            k2 = IA(ii+1) - 1
            do jj = k1,k2
                if (JA(jj) <= ii) nnz = nnz + 1
            end do
        end do

        allocate(L_AA(nnz))
        allocate(L_JA(nnz))
        allocate(L_IA(n+1))
        allocate(L_D0(n))

        kk = 0
        do ii = 1, n
            L_IA(ii) = kk + 1
            k1 = IA(ii)
            k2 = IA(ii+1) - 1
            do jj = k1, k2
                if (JA(jj) <= ii) then
                    kk = kk + 1
                    L_AA(kk) = AA(jj)
                    L_JA(kk) = JA(jj)
                end if
            end do
        end do
        L_IA(n+1) = nnz + 1  ! last element points past the end

        do ii = 1, n
            k1 = L_IA(ii)
            k2 = L_IA(ii+1)-1
            do jj = k1, k2
                jcol = L_JA(jj)
                s = L_AA(jj)
                do kk1 = k1, jj-1
                    kcol = L_JA(kk1) ! This is k
                    ! We need L(j,k), which is L(jcol, kcol). Search row jcol for column kcol.
                    do kk2 = L_IA(jcol), L_IA(jcol+1)-1
                        if (L_JA(kk2) == kcol) then
                            s = s - L_AA(kk1) * L_AA(kk2) ! L(i,k) * L(j,k)
                            exit
                        end if
                    end do
                end do

                if (ii == jcol) then
                    s = max(s, 1.0d-20)
                    L_D0(ii) = sqrt(s)
                    L_AA(jj) = L_D0(ii)
                else
                    ! This is L(i,j). We need L(j,j) which is L_D0(jcol)
                    L_AA(jj) = s / L_D0(jcol)
                end if
            end do
        end do
    end subroutine

    !-----------------------------------------------------------------
    ! Performs ILU(0) factorization of sparse matrix A in CSR format:
    ! A ≈ L * U   with unit diagonal in L
    ! Pattern of L and U follows pattern of A
    !-----------------------------------------------------------------
    subroutine ILU0_CSR(AA, JA, IA, L_AA, L_JA, L_IA, U_AA, U_JA, U_IA)
        implicit none
        real(8), intent(in) :: AA(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8), allocatable, intent(out) :: L_AA(:), U_AA(:)
        integer, allocatable, intent(out) :: L_JA(:), L_IA(:), U_JA(:), U_IA(:)

        integer :: n, i, j, k, p, p2
        real(8), allocatable :: w(:)
        integer, allocatable :: iw(:)

        n = size(IA) - 1
        allocate(L_IA(n+1), U_IA(n+1))
        L_IA = IA; U_IA = IA
        allocate(L_AA(size(AA)), L_JA(size(JA)), U_AA(size(AA)), U_JA(size(JA)))
        L_AA = AA; L_JA = JA; U_AA = AA; U_JA = JA

        allocate(w(n), iw(n))
        iw = 0

        do i = 1, n
            ! Scatter row i into work vector w
            do p = IA(i), IA(i+1)-1
                w(JA(p)) = AA(p)
                iw(JA(p)) = i
            end do

            ! Update row i based on previous rows
            do p = IA(i), IA(i+1)-1
                k = JA(p)
                if (k >= i) exit ! Assumes sorted columns
                
                ! L(i,k) = A(i,k) / U(k,k)
                w(k) = w(k) / U_AA(U_IA(k+1)-1) ! U(k,k) is the last element in row k of U

                ! Update remaining elements in row i
                do p2 = U_IA(k), U_IA(k+1)-1
                    j = U_JA(p2)
                    if (iw(j) == i) then ! If A(i,j) is non-zero
                        w(j) = w(j) - w(k) * U_AA(p2)
                    end if
                end do
            end do

            ! Gather updated row back into L and U
            do p = IA(i), IA(i+1)-1
                j = JA(p)
                if (j < i) then
                    L_AA(p) = w(j)
                elseif (j == i) then
                    L_AA(p) = 1.0_dp
                    U_AA(p) = max(w(j), 1.0d-20)
                else ! j > i
                    U_AA(p) = w(j)
                end if
            end do

            ! Reset work vector's sparsity tracker
            do p = IA(i), IA(i+1)-1
                iw(JA(p)) = 0
            end do
        end do

        deallocate(w, iw)
    end subroutine

    !-------------------------
    ! Function to compute diagonal on input CSR matrix A
    ! Returns diag, the diagonal elements of A
    !-------------------------
    subroutine Jacobi_CSR(AA, JA, IA, diag)
        implicit none
        ! Inputs
        real(8), intent(in) :: AA(:)
        integer, intent(in) :: JA(:), IA(:)
        ! Output
        real(8), intent(out), allocatable :: diag(:)

        integer :: n, ii, jj

        n = size(IA) - 1
        allocate(diag(n))

        do ii = 1, n
            ! Loop over row ii
            do jj = IA(ii), IA(ii+1)-1
                if (JA(jj) == ii) then
                    diag(ii) = AA(jj)
                    exit
                end if
            end do
        end do
    end subroutine

    !-------------------------
    ! Function to compute d = M^-1 r using preconditioning matrix M from IC(0) factorization
    ! L_AA, L_JA, L_IA (CSR format) define preconditioning matrix M computed by Cholesky Decomposition
    ! Returns search vector d
    !-------------------------
    function PCG_Cholesky_CSR(L_AA, L_JA, L_IA, r) result(d)
        real(8), intent(in) :: L_AA(:), r(:)
        integer, intent(in) :: L_JA(:), L_IA(:)
        real(8), allocatable :: d(:), y(:)

        !Incomplete  Cholesky preconditioning: solve M d = r  with M = L*L^T
        allocate(y(size(L_IA)-1))
        call forward_solver_CSR(L_AA,L_JA,L_IA,r,y)
        allocate(d(size(L_IA)-1))
        call backward_solver_CSR(L_AA,L_JA,L_IA,y,d)
    end function

    !-------------------------
    ! Returns search vector d
    !-------------------------
    function PCG_ILU_CSR(L_AA, L_JA, L_IA, U_AA, U_JA, U_IA, r) result(d)
        real(8), intent(in) :: L_AA(:), U_AA(:), r(:)
        integer, intent(in) :: L_JA(:), L_IA(:), U_JA(:), U_IA(:)
        real(8), allocatable :: d(:), y(:)
        
        !Incomplete LU preconditioning: solve M d = r  with M = L*U
        allocate(y(size(L_IA)-1))
        call forward_solver_CSR(L_AA,L_JA,L_IA,r,y)
        allocate(d(size(L_IA)-1))
        call backward_solver_upper_CSR(U_AA,U_JA,U_IA,y,d)

    end function

    !-------------------------
    ! Forward solve: L y = r  (L stored in CSR; diagonal is last element per row)
    !-------------------------
    subroutine forward_solver_CSR(L_AA, L_JA, L_IA, r0, y0)
        implicit none
        real(8), intent(in) :: L_AA(:), r0(:)
        integer, intent(in) :: L_JA(:), L_IA(:)
        real(8), intent(out) :: y0(:)
        integer :: ii, k1, k2, jj
        real(8) :: diag_val, sum_val

        ! Solve L*y = r where L is lower triangular stored in CSR
        ! L has the structure from Cholesky decomposition: diagonal is stored inline
        
        do ii = 1, size(L_IA) - 1
            k1 = L_IA(ii)
            k2 = L_IA(ii+1) - 1
            
            ! Find diagonal and compute forward substitution
            diag_val = 0.0d0
            sum_val = 0.0d0
            
            do jj = k1, k2
                if (L_JA(jj) == ii) then
                    diag_val = L_AA(jj)
                else if (L_JA(jj) < ii) then
                    sum_val = sum_val + L_AA(jj) * y0(L_JA(jj))
                end if
            end do
            
            if (abs(diag_val) < 1.0d-20) then
                print *, "ERROR: Near-zero diagonal in forward_solver at row", ii, ":", diag_val
                y0(ii) = 0.0d0
            else
                y0(ii) = (r0(ii) - sum_val) / diag_val
            end if
        end do
        return
    end subroutine

    !-------------------------
    ! Backward solve: L^T d = y  (L stored in CSR; diagonal is stored inline)
    !-------------------------
    subroutine backward_solver_CSR(L_AA,L_JA,L_IA,y0,d0)
        implicit none
        real(8), intent(in) :: L_AA(:), y0(:)
        integer, intent(in) :: L_JA(:), L_IA(:)
        real(8),intent(out) :: d0(:)

        integer :: i, j, k, n
        real(8) :: diag_val

        ! Solve L^T * d = y where L is stored in CSR format
        
        n = size(L_IA) - 1
        d0 = y0

        do i = n, 1, -1
            ! Find diagonal L(i,i)
            diag_val = 0.0d0
            do k = L_IA(i), L_IA(i+1)-1
                if (L_JA(k) == i) then
                    diag_val = L_AA(k)
                    exit
                end if
            end do
            
            if (abs(diag_val) < 1.0d-20) then
                d0(i) = 0.0d0
            else
                d0(i) = d0(i) / diag_val
            end if

            ! Update d0(j) for j < i using L(i,j)
            ! This corresponds to subtracting L(i,j)*d(i) from row j of the system
            do k = L_IA(i), L_IA(i+1)-1
                j = L_JA(k)
                if (j < i) then
                    d0(j) = d0(j) - L_AA(k) * d0(i)
                end if
            end do
        end do
        
        return
    end subroutine

    !-------------------------
    ! Backward solve for upper triangular: U*d = y
    ! U is upper triangular stored in CSR; solve by backward substitution
    !-------------------------
    subroutine backward_solver_upper_CSR(U_AA, U_JA, U_IA, y0, d0)
        implicit none
        real(8), intent(in) :: U_AA(:), y0(:)
        integer, intent(in) :: U_JA(:), U_IA(:)
        real(8), intent(out) :: d0(size(U_IA)-1)

        integer :: ii, jj, k1, k2, n
        real(8) :: diag_val, sum_val

        ! Solve U*d = y where U is upper triangular stored in CSR
        ! U(i,j) with j >= i are stored
        
        n = size(U_IA) - 1
        d0 = y0

        ! Backward substitution: start from last row and go backwards
        do ii = n, 1, -1
            k1 = U_IA(ii)
            k2 = U_IA(ii+1) - 1
            
            ! Find diagonal element U(ii,ii)
            diag_val = 0.0d0
            sum_val = 0.0d0
            
            do jj = k1, k2
                if (U_JA(jj) == ii) then
                    diag_val = U_AA(jj)
                else if (U_JA(jj) > ii) then
                    ! Upper entries: U(ii, jj) where jj > ii
                    sum_val = sum_val + U_AA(jj) * d0(U_JA(jj))
                end if
            end do
            
            if (abs(diag_val) < 1.0d-20) then
                print *, "ERROR: Near-zero diagonal in backward_solver_upper at row", ii, ":", diag_val
                d0(ii) = 0.0d0
            else
                d0(ii) = (d0(ii) - sum_val) / diag_val
            end if
        end do

        return
    end subroutine

end module