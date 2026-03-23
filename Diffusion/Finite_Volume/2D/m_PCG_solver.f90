
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
! 11/11/25     T. Charlton    Implemented 2D Application                -!
!------------------------------------------------------------------------! 

module m_PCG_solver
implicit none
    !-----------------------------
    ! Enumerated preconditioner types
    !-----------------------------
    integer, parameter :: PCG_PRECON_NONE = 0
    integer, parameter :: PCG_PRECON_CHOLESKY = 1
    integer, parameter :: PCG_PRECON_ILU = 2
    integer, parameter :: PCG_PRECON_JACOBI = 3

contains 

    !-------------------------
    ! Principle Conjugation Gradient (PCG) Algorithm driver, solving Ax=b
    !-------------------------
    subroutine PCG_algorithm(AA, JA, IA, x, b, PCG_mode, max_iter)
        implicit none 

        real(8) :: alpha, beta
        integer :: ii
        real(8), intent(in) :: AA(:), b(:)
        integer, intent(in) :: JA(:), IA(:)
        integer, intent(in) :: PCG_mode, max_iter
        real(8), intent(out) :: x(:)

        real(8) :: rho_old, rho_new, denom !residual_norm
        real(8), allocatable :: L_AA(:), U_AA(:), diag(:), r(:), d(:), z(:), q(:)
        integer, allocatable :: L_IA(:), L_JA(:), U_IA(:), U_JA(:)
        logical :: nan_detected

        allocate(r(size(b)), d(size(b)), z(size(b)), q(size(b)))

        r = b - CSR_dot_product(AA, JA, IA, x)

        ! apply preconditioner to get initial z = M^-1 r
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

            x = x + alpha*d
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

            if (abs(rho_old) < 1.0d-21) then
                !print *, "WARNING: rho_old too small in iteration", ii, ":", rho_old
                exit
            end if
            
            !print*, ii, r(1), r(50), r(2000), r(4000)
            
            beta = rho_new / rho_old
            d = z + beta*d
            rho_old = rho_new
        end do

        if (any(isnan(x))) then
            print *, "ERROR: NaN detected in final solution x"
        end if

    end subroutine

    !-------------------------
    ! CSR Matvec multiplication function
    !-------------------------
    function CSR_dot_product(AA,JA,IA,x) result(Y)
        implicit none
        real(8), intent(in) :: AA(:), x(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8) :: Y(SIZE(x))
        integer :: ii, k1, k2

        do ii = 1,SIZE(x)
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