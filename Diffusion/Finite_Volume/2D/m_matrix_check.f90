module m_matrix_check
  implicit none
  private
  public :: check_matrix_SPD_CSR, check_matrix_symmetry_CSR, check_matrix_positive_definite_CSR

contains

  !---------------------------------------------------------------
  ! Check if a CSR matrix is Symmetric Positive Definite (SPD)
  !---------------------------------------------------------------
  subroutine check_matrix_SPD_CSR(ia, ja, a, is_spd, is_symmetric, is_positive_definite, verbose)
    implicit none
    integer, intent(in) :: ia(:), ja(:)
    real(8), intent(in) :: a(:)
    logical, intent(out) :: is_spd, is_symmetric, is_positive_definite
    logical, intent(in), optional :: verbose

    logical :: verb
    integer :: n

    verb = .false.
    if (present(verbose)) verb = verbose

    n = size(ia) - 1

    ! Check symmetry first
    call check_matrix_symmetry_CSR(ia, ja, a, is_symmetric, verb)

    if (.not. is_symmetric) then
      is_spd = .false.
      is_positive_definite = .false.
      if (verb) print *, "Matrix is not symmetric. Cannot be SPD."
      return
    end if

    ! Check positive definiteness
    call check_matrix_positive_definite_CSR(ia, ja, a, is_positive_definite, verb)

    is_spd = is_symmetric .and. is_positive_definite

    if (verb) then
      if (is_spd) then
        print *, "Matrix IS Symmetric Positive Definite (SPD)"
      else
        print *, "Matrix is NOT Symmetric Positive Definite (SPD)"
      end if
    end if

  end subroutine check_matrix_SPD_CSR

  !---------------------------------------------------------------
  ! Check if a CSR matrix is symmetric
  !---------------------------------------------------------------
  subroutine check_matrix_symmetry_CSR(ia, ja, a, is_symmetric, verbose)
    implicit none
    integer, intent(in) :: ia(:), ja(:)
    real(8), intent(in) :: a(:)
    logical, intent(out) :: is_symmetric
    logical, intent(in), optional :: verbose

    logical :: verb
    integer :: n, i, j, k, k_j
    real(8) :: a_ij, a_ji, tol
    logical :: found

    verb = .false.
    if (present(verbose)) verb = verbose

    tol = 1.0d-10
    n = size(ia) - 1
    is_symmetric = .true.

    do i = 1, n
      do k = ia(i), ia(i+1) - 1
        j = ja(k)
        a_ij = a(k)

        ! Find a(j,i)
        found = .false.
        do k_j = ia(j), ia(j+1) - 1
          if (ja(k_j) == i) then
            a_ji = a(k_j)
            found = .true.
            exit
          end if
        end do

        if (.not. found) then
          if (verb) print *, "Symmetry violation: a(", i, ",", j, ") = ", a_ij, " but a(", j, ",", i, ") missing"
          is_symmetric = .false.
          return
        end if

        if (abs(a_ij - a_ji) > tol) then
          if (verb) print *, "Symmetry violation: a(", i, ",", j, ") = ", a_ij, " but a(", j, ",", i, ") = ", a_ji
          is_symmetric = .false.
          return
        end if
      end do
    end do

    ! if (verb) then
    !   if (is_symmetric) then
    !     print *, "Matrix is SYMMETRIC"
    !   else
    !     print *, "Matrix is NOT SYMMETRIC"
    !   end if
    ! end if

  end subroutine check_matrix_symmetry_CSR

  !---------------------------------------------------------------
  ! Check if a symmetric CSR matrix is positive definite
  ! Uses Cholesky factorization attempt (Incomplete Cholesky IC(0))
  !---------------------------------------------------------------
  subroutine check_matrix_positive_definite_CSR(ia, ja, a, is_positive_definite, verbose)
    implicit none
    integer, intent(in) :: ia(:), ja(:)
    real(8), intent(in) :: a(:)
    logical, intent(out) :: is_positive_definite
    logical, intent(in), optional :: verbose

    logical :: verb
    real(8), allocatable :: L_aa(:), L_D(:)
    integer, allocatable :: L_ja(:), L_ia(:)
    real(8) :: s, min_diag, max_diag
    integer :: ii, jj, kk, k1, k2, nnz, jcol, kcol, kk1, kk2, n
    logical :: breakdown

    verb = .false.
    if (present(verbose)) verb = verbose

    n = size(ia) - 1
    is_positive_definite = .true.
    breakdown = .false.

    ! Count nonzeros in lower triangle
    nnz = 0
    do ii = 1, n
      k1 = ia(ii)
      k2 = ia(ii+1) - 1
      do jj = k1, k2
        if (ja(jj) <= ii) nnz = nnz + 1
      end do
    end do

    allocate(L_aa(nnz))
    allocate(L_ja(nnz))
    allocate(L_ia(n+1))
    allocate(L_D(n))

    ! Extract lower triangle
    kk = 0
    do ii = 1, n
      L_ia(ii) = kk + 1
      k1 = ia(ii)
      k2 = ia(ii+1) - 1
      do jj = k1, k2
        if (ja(jj) <= ii) then
          kk = kk + 1
          L_aa(kk) = a(jj)
          L_ja(kk) = ja(jj)
        end if
      end do
    end do
    L_ia(n+1) = nnz + 1

    min_diag = huge(1.0d0)
    max_diag = 0.0d0

    ! Attempt Incomplete Cholesky factorization
    do ii = 1, n
      k1 = L_ia(ii)
      k2 = L_ia(ii+1) - 1
      do jj = k1, k2
        jcol = L_ja(jj)
        s = L_aa(jj)
        do kk1 = k1, jj-1
          kcol = L_ja(kk1)
          do kk2 = L_ia(kcol), L_ia(kcol+1)-1
            if (L_ja(kk2) == jcol) then
              s = s - L_aa(kk1)*L_aa(kk2)
              exit
            end if
          end do
        end do

        if (ii == jcol) then
          ! Diagonal element
          if (s <= 0.0d0) then
            if (verb) print *, "Non-positive diagonal at row ", ii, ": s = ", s
            breakdown = .true.
            is_positive_definite = .false.
            exit
          end if
          L_D(ii) = sqrt(s)
          min_diag = min(min_diag, L_D(ii))
          max_diag = max(max_diag, L_D(ii))
        else
          ! Off-diagonal element
          if (abs(L_D(jcol)) < 1.0d-20) then
            if (verb) print *, "Near-zero diagonal at column ", jcol, " encountered"
            breakdown = .true.
            is_positive_definite = .false.
            exit
          end if
        end if
      end do

      if (breakdown) exit
    end do

    if (verb) then
      if (is_positive_definite) then
        ! print *, "Matrix is POSITIVE DEFINITE"
        print *, "  Diagonal range: [", min_diag, ", ", max_diag, "]"
        print *, "  Condition number estimate: ", max_diag / min_diag
      else
        print *, "Matrix is NOT POSITIVE DEFINITE (Cholesky factorization failed)"
      end if
    end if

    deallocate(L_aa, L_ja, L_ia, L_D)

  end subroutine check_matrix_positive_definite_CSR

  !---------------------------------------------------------------
  ! Compute diagonal of CSR matrix
  !---------------------------------------------------------------
  subroutine get_matrix_diagonal_CSR(ia, ja, a, diag)
    implicit none
    integer, intent(in) :: ia(:), ja(:)
    real(8), intent(in) :: a(:)
    real(8), intent(out) :: diag(:)

    integer :: i, k

    diag = 0.0d0

    do i = 1, size(ia) - 1
      do k = ia(i), ia(i+1) - 1
        if (ja(k) == i) then
          diag(i) = a(k)
          exit
        end if
      end do
    end do

  end subroutine get_matrix_diagonal_CSR

  !---------------------------------------------------------------
  ! Check if all diagonal elements are positive
  !---------------------------------------------------------------
  subroutine check_diagonal_positive_CSR(ia, ja, a, all_positive, min_diag, max_diag, verbose)
    implicit none
    integer, intent(in) :: ia(:), ja(:)
    real(8), intent(in) :: a(:)
    logical, intent(out) :: all_positive
    real(8), intent(out) :: min_diag, max_diag
    logical, intent(in), optional :: verbose

    logical :: verb
    real(8), allocatable :: diag(:)
    integer :: i, n

    verb = .false.
    if (present(verbose)) verb = verbose

    n = size(ia) - 1
    allocate(diag(n))

    call get_matrix_diagonal_CSR(ia, ja, a, diag)

    min_diag = minval(diag)
    max_diag = maxval(diag)
    all_positive = all(diag > 0.0d0)

    if (verb) then
      print *, "Diagonal element range: [", min_diag, ", ", max_diag, "]"
      if (all_positive) then
        print *, "All diagonal elements are POSITIVE"
      else
        print *, "Some diagonal elements are NON-POSITIVE"
        do i = 1, n
          if (diag(i) <= 0.0d0) print *, "  Row ", i, ": diag = ", diag(i)
        end do
      end if
    end if

    deallocate(diag)

  end subroutine check_diagonal_positive_CSR

end module m_matrix_check
