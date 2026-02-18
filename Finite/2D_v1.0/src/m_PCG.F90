
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
use m_finite_elements
use m_material
use m_quadrature
use m_types
implicit none
    integer, parameter :: PCG_PRECON_NONE = 0
    integer, parameter :: PCG_PRECON_CHOLESKY = 1
    integer, parameter :: PCG_PRECON_ILU = 2
    integer, parameter :: PCG_PRECON_JACOBI = 3
    integer, parameter :: PETsc = 4

contains 

    subroutine assemble_PCG_matrix(MAT_DATA, A_MAT_PCG, mesh, FE, Quad, NGRP, MATIDS)
        type(t_material),       intent(in)  :: MAT_DATA(:)
        type(t_mesh),           intent(in)  :: mesh
        type(t_finite),         intent(in)  :: FE
        type(t_quadrature),     intent(in)  :: Quad
        integer,                intent(in)  :: NGRP
        integer,                intent(in)  :: MATIDS(:)
        type(t_mat), allocatable,  intent(out) :: A_MAT_PCG(:)
        
        integer :: ee, g, i, j, q, mat_id, row, col
        integer, allocatable :: nnz(:)
        real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV, val

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
            call PCG_MAT_INIT(A_MAT_PCG, mesh, nnz, g)
            A_MAT_PCG(g)%row_ptr(1) = 1
            do i = 1, mesh%n_nodes
                A_MAT_PCG(g)%row_ptr(i+1) = A_MAT_PCG(g)%row_ptr(i) + nnz(i)
            end do
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
                        row = mesh%elems(ee, i)
                        do j = 1, FE%n_basis
                            col = mesh%elems(ee, j)
                            
                            val = (MAT_DATA(mat_id)%D(g) * (dN_dx(i)*dN_dx(j) + dN_dy(i)*dN_dy(j)) + &
                                MAT_DATA(mat_id)%SigmaR(g) * FE%N_mat(q,i,j)) * dV
                                
                            call PCG_MAT_ALLOCATION(A_MAT_PCG(g), row, col, val)
                        end do
                    end do
                end do
            end do
        end do

        do g = 1, NGRP
            call PCG_COMPRESS_CSR(A_MAT_PCG(g), mesh%n_nodes)
        end do

        deallocate(nnz)
    end subroutine assemble_PCG_matrix

    subroutine assemble_PCG_source_vec(B_VEC, mesh, FE, Quad, NGRP, S_ext)
        type(t_mesh),           intent(in)  :: mesh
        type(t_finite),         intent(in)  :: FE
        type(t_quadrature),     intent(in)  :: Quad
        integer,                intent(in)  :: NGRP
        real(dp),               intent(in)  :: S_ext(:,:) 
        type(t_VEC), allocatable,  intent(inout) :: B_VEC(:)

        integer  :: ee, g, i, q, row
        real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV, local_val

        if (allocated(B_VEC)) deallocate(B_VEC)
        allocate(B_VEC(NGRP))

        do g = 1, NGRP
            call PCG_VEC_INIT(mesh%n_nodes, B_VEC(g))
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
                        row = mesh%elems(ee, i)
                        local_val = S_ext(ee, g) * FE%N(q, i) * dV
                        call PCG_VEC_ALLOCATION(B_VEC(g), row, local_val)
                    end do
                end do
            end do
        end do
    end subroutine assemble_PCG_source_vec

    subroutine PCG_MAT_INIT(A_MAT_PCG, mesh, nnz_array, g)
        type(t_mesh), intent(in) :: mesh
        type(t_mat),  intent(inout) :: A_MAT_PCG(:)
        integer,      intent(in)    :: nnz_array(:)
        integer,      intent(in)    :: g
        integer :: total_nnz

        total_nnz = sum(nnz_array)

        if (allocated(A_MAT_PCG(g)%val)) deallocate(A_MAT_PCG(g)%val)
        if (allocated(A_MAT_PCG(g)%col)) deallocate(A_MAT_PCG(g)%col)
        if (allocated(A_MAT_PCG(g)%row_ptr)) deallocate(A_MAT_PCG(g)%row_ptr)

        allocate(A_MAT_PCG(g)%val(total_nnz))
        allocate(A_MAT_PCG(g)%col(total_nnz))
        allocate(A_MAT_PCG(g)%row_ptr(mesh%n_nodes + 1))
        
        A_MAT_PCG(g)%val = 0.0_dp
        A_MAT_PCG(g)%col = 0
        A_MAT_PCG(g)%row_ptr = 0
    end subroutine PCG_MAT_INIT

    subroutine PCG_MAT_ALLOCATION(A_MAT_PCG, row, col, val)
        type(t_mat), intent(inout) :: A_MAT_PCG
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
        type(t_mat), intent(inout) :: A_MAT_PCG
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

    subroutine assemble_multigroup_fission_pcg(Q_PCG, mesh, FE, Quad, materials, X_PCG, n_groups, k_eff)
        type(t_vec), intent(inout)      :: Q_PCG
        type(t_mesh), intent(in)        :: mesh
        type(t_finite), intent(in)      :: FE
        type(t_quadrature), intent(in)  :: Quad
        type(t_material), intent(in)    :: materials(:)
        type(t_vec), intent(in)         :: X_PCG(:)
        integer, intent(in)             :: n_groups
        real(dp), intent(in)            :: k_eff

        integer :: ee, g, i, q, mat_id
        real(dp) :: detJ, dV, local_phi_fission, inv_k
        real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), phi_nodes(FE%n_basis)
        integer  :: idx(FE%n_basis)
        
        inv_k = 1.0_dp / k_eff
        Q_PCG%vec = 0.0_dp

        do ee = 1, mesh%n_elems
            mat_id = mesh%mats(ee)
            idx = mesh%elems(ee, :)
            do i = 1, FE%n_basis; elem_coords(i, :) = mesh%nodes(idx(i), :); end do

            do q = 1, Quad%NoPoints
                call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
                dV = detJ * Quad%W(q)
                
                local_phi_fission = 0.0_dp
                do g = 1, n_groups
                    phi_nodes = X_PCG(g)%vec(idx)
                    local_phi_fission = local_phi_fission + materials(mat_id)%nuSigF(g) * dot_product(FE%N(q, :), phi_nodes)
                end do

                do i = 1, FE%n_basis
                    Q_PCG%vec(idx(i)) = Q_PCG%vec(idx(i)) + local_phi_fission * FE%N(q, i) * dV
                end do
            end do
        end do
    end subroutine

    subroutine calculate_total_production_pcg(total_prod, X_PCG, mesh, FE, Quad, materials, n_groups)
        real(dp), intent(out)           :: total_prod
        type(t_vec), intent(in)         :: X_PCG(:)
        type(t_mesh), intent(in)        :: mesh
        type(t_finite), intent(in)      :: FE
        type(t_quadrature), intent(in)  :: Quad
        type(t_material), intent(in)    :: materials(:)
        integer, intent(in)             :: n_groups

        integer  :: ee, g, i, q, mat_id
        integer  :: idx(FE%n_basis)
        real(dp) :: detJ, dV, local_phi, elem_phi(FE%n_basis)
        real(dp) :: elem_coords(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis)

        total_prod = 0.0_dp

        do ee = 1, mesh%n_elems
            mat_id = mesh%mats(ee)
            idx = mesh%elems(ee, :)
            
            do i = 1, FE%n_basis
                elem_coords(i, :) = mesh%nodes(idx(i), :)
            end do

            do q = 1, Quad%NoPoints
                call GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
                dV = detJ * Quad%W(q)
                
                do g = 1, n_groups
                    elem_phi = X_PCG(g)%vec(idx)
                    
                    local_phi = dot_product(FE%N(q, :), elem_phi)
                    
                    total_prod = total_prod + materials(mat_id)%nuSigF(g) * local_phi * dV
                end do
            end do
        end do
    end subroutine calculate_total_production_pcg

    !-------------------------
    ! Principle Conjugation Gradient (PCG) Algorithm driver, solving Ax=b
    !-------------------------
    subroutine CGSolve(A_MAT_PCG, x, b, PCG_mode, max_iter)
        implicit none 
        type(t_mat), intent(in) :: A_MAT_PCG
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
            case(PCG_PRECON_NONE)
                z = r
            case(PCG_PRECON_CHOLESKY)
                call Cholesky_CSR(AA, JA, IA, L_AA, L_JA, L_IA)
                z = PCG_Cholesky_CSR(L_AA, L_JA, L_IA, r)
            case(PCG_PRECON_ILU)
                call ILU0_CSR(AA, JA, IA, L_AA, L_JA, L_IA, U_AA, U_JA, U_IA)
                z = PCG_ILU_CSR(L_AA, L_JA, L_IA, U_AA, U_JA, U_IA, r)
            case(PCG_PRECON_JACOBI)
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
                case(PCG_PRECON_NONE)
                    z = r
                case(PCG_PRECON_CHOLESKY)
                    z = PCG_Cholesky_CSR(L_AA, L_JA, L_IA, r)
                case(PCG_PRECON_ILU)
                    z = PCG_ILU_CSR(L_AA, L_JA, L_IA, U_AA, U_JA, U_IA, r)
                case(PCG_PRECON_JACOBI)
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
                    kcol = L_JA(kk1)
                    do kk2 = L_IA(kcol), L_IA(kcol+1)-1
                        if (L_JA(kk2) == jcol) then
                            s = s - L_AA(kk1)*L_AA(kk2)
                            exit
                        end if
                    end do
                end do

                if (ii == jcol) then
                    ! diagonal
                    s = max(s, 1.0d-20)
                    L_D0(ii) = sqrt(s)
                    L_AA(jj) = L_D0(ii)
                else
                    ! off-diagonal
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

        integer :: n, nnzL, nnzU
        integer :: ii, jj, k1, k2, kk
        integer :: col
        real(8) :: s

        n = size(IA) - 1

        !--------------------------------------------------
        ! Count number of nonzeros for L and U
        ! L: lower-triangular + diagonal
        ! U: upper-triangular (including diagonal)
        !--------------------------------------------------
        nnzL = 0
        nnzU = 0
        do ii = 1, n
            do jj = IA(ii), IA(ii+1)-1
                if (JA(jj) <= ii) then
                    nnzL = nnzL + 1
                end if
                if (JA(jj) >= ii) then
                    nnzU = nnzU + 1
                end if
            end do
        end do

        allocate(L_AA(nnzL), L_JA(nnzL), L_IA(n+1))
        allocate(U_AA(nnzU), U_JA(nnzU), U_IA(n+1))

        kk = 0
        do ii = 1, n
            L_IA(ii) = kk + 1
            do jj = IA(ii), IA(ii+1)-1
                if (JA(jj) <= ii) then
                    kk = kk + 1
                    L_JA(kk) = JA(jj)
                    L_AA(kk) = AA(jj)
                end if
            end do
        end do
        L_IA(n+1) = kk + 1

        kk = 0
        do ii = 1, n
            U_IA(ii) = kk + 1
            do jj = IA(ii), IA(ii+1)-1
                if (JA(jj) >= ii) then
                    kk = kk + 1
                    U_JA(kk) = JA(jj)
                    U_AA(kk) = AA(jj)
                end if
            end do
        end do
        U_IA(n+1) = kk + 1

        do ii = 1, n
            ! Process L row
            do jj = L_IA(ii), L_IA(ii+1)-1
                col = L_JA(jj)
                s = L_AA(jj)
                ! Subtract previous contributions
                do k1 = L_IA(ii), jj-1
                    if (L_JA(k1) < col) cycle
                    ! Find matching U element
                    do k2 = U_IA(L_JA(k1)), U_IA(L_JA(k1)+1)-1
                        if (U_JA(k2) == col) then
                            s = s - L_AA(k1)*U_AA(k2)
                            exit
                        end if
                    end do
                end do
                ! Divide by diagonal of U
                if (col /= ii) then
                    ! Find U diagonal
                    do k2 = U_IA(col), U_IA(col+1)-1
                        if (U_JA(k2) == col) then
                            L_AA(jj) = s / U_AA(k2)
                            exit
                        end if
                    end do
                else
                    L_AA(jj) = s
                end if
            end do

            ! Process U row
            do jj = U_IA(ii), U_IA(ii+1)-1
                col = U_JA(jj)
                s = U_AA(jj)
                do k1 = L_IA(ii), L_IA(ii+1)-1
                    if (L_JA(k1) < ii) then
                        do k2 = U_IA(L_JA(k1)), U_IA(L_JA(k1)+1)-1
                            if (U_JA(k2) == col) then
                                s = s - L_AA(k1)*U_AA(k2)
                                exit
                            end if
                        end do
                    end if
                end do
                if (col == ii) s = max(s, 1.0d-20)
                U_AA(jj) = s
            end do
        end do
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
        real(8),intent(out) :: d0(size(L_IA)-1)

        integer :: ii, jj, k1, k2, n
        real(8) :: diag_val

        ! Solve L^T * d = y where L is stored in CSR format
        ! For the transpose: L^T has L(i,j) values at position (j,i)
        ! We need to iterate backwards and handle the transpose implicitly
        
        n = size(L_IA) - 1
        d0 = y0

        ! Backward substitution for L^T
        do ii = n, 1, -1
            k1 = L_IA(ii)
            k2 = L_IA(ii+1) - 1
            
            ! Find diagonal element in row ii
            diag_val = 0.0d0
            do jj = k1, k2
                if (L_JA(jj) == ii) then
                    diag_val = L_AA(jj)
                    exit
                end if
            end do
            
            if (abs(diag_val) < 1.0d-20) then
                print *, "ERROR: Near-zero diagonal in backward_solver at row", ii, ":", diag_val
                d0(ii) = 0.0d0
            else
                d0(ii) = d0(ii) / diag_val
            end if

            ! Subtract contributions from L^T(j,i) where j > i
            ! These correspond to L(i,j) entries where j > i
            ! We need to update d0(jj) for jj > ii using values L(ii,jj)
            do jj = k1, k2
                if (L_JA(jj) > ii) then
                    ! This entry is L(ii, L_JA(jj)), and in L^T it's at (L_JA(jj), ii)
                    ! But since we're doing backward substitution, we need forward entries
                    ! This approach won't work - we need to iterate over future rows
                end if
            end do
        end do

        ! Better approach: iterate over all rows to find L^T contributions
        do ii = n, 1, -1
            if (abs(d0(ii)) < 1.0d-20) cycle
            
            ! Find diagonal in row ii
            k1 = L_IA(ii)
            k2 = L_IA(ii+1) - 1
            diag_val = 0.0d0
            
            do jj = k1, k2
                if (L_JA(jj) == ii) then
                    diag_val = L_AA(jj)
                    exit
                end if
            end do
            
            if (abs(diag_val) < 1.0d-20) then
                if (ii > 1) then  ! Not catastrophic for last rows
                    print *, "WARNING: Near-zero diagonal in backward_solver at row", ii
                end if
                cycle
            end if
            
            ! Update d0(ii)
            d0(ii) = d0(ii) / diag_val
            
            ! Update rows jj > ii using the L(ii,jj) entries (which become L^T(jj,ii))
            do jj = k1, k2
                if (L_JA(jj) > ii) then
                    d0(L_JA(jj)) = d0(L_JA(jj)) - L_AA(jj) * d0(ii)
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