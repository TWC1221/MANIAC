module m_diffusion_matrix
  use m_constants
  implicit none
  contains

  subroutine assemble_2D_diffusion_matrix(nx, ny, x_domain, y_domain, D, Sigma_a, ia, ja, a, dx, dy)
    implicit none
    !---------------------------------------------------------------
    ! Assembles CSR matrix for 2D diffusion operator:
    !   -D*(d2/dx2 + d2/dy2) + Sigma_a
    ! Dirichlet boundary conditions (phi=0) on all edges.
    !
    ! Input:
    !   nx, ny     - grid points (including boundaries)
    !   dx, dy     - grid spacing
    !   D          - diffusion coefficient
    !   Sigma_a    - absorption coefficient
    !
    ! Output:
    !   ia(:)      - kk pointer (CSR)
    !   ja(:)      - column indices (CSR)
    !   a(:)       - nonzero values (CSR)
    !
    ! Notes:
    !   Lexicographic ordering: k = ii + (jj-1)*nx
    !---------------------------------------------------------------
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: x_domain, y_domain, D, Sigma_a
    real(8), intent(out) :: dx, dy
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:)

    integer :: N, ii, jj, kk, ptr, est_nnz
    real(8) :: aL, aR, aB, aT, aC

    dx = x_domain / real(nx) ; dy = y_domain / real(ny)

    N = nx * ny
    est_nnz = 5*N

    allocate(ia(N+1), ja(est_nnz), a(est_nnz))
    ptr = 1
    ia(1) = 1

    do jj = 1, ny
      do ii = 1, nx
        kk = ii + (jj-1)*nx

        ! Base interior coefficients
        aL = -D * dy / dx
        aR = -D * dy / dx
        aB = -D * dx / dy
        aT = -D * dx / dy

        ! Multiply absorption by cell area (finite-volume form)
        aC = Sigma_a*dx*dy - (aL + aR + aB + aT)

        !-----------------------------------------------------
        ! Apply Dirichlet BCs via zero-value ghost cells
        ! Remove outside neighbors and add to diagonal
        !-----------------------------------------------------

        ! Left boundary face
        if (ii == 1) then
            aC = aC - aL
            aL = 0.0
        end if

        ! Right boundary
        if (ii == nx) then
            aC = aC - aR
            aR = 0.0
        end if

        ! Bottom
        if (jj == 1) then
            aC = aC - aB
            aB = 0.0
        end if

        ! Top
        if (jj == ny) then
            aC = aC - aT
            aT = 0.0
        end if

        !-----------------------------------------------------
        ! Write row in CSR
        !-----------------------------------------------------

        ! Center
        ja(ptr) = kk
        a(ptr)  = aC
        ptr = ptr + 1

        ! Left neighbor
        if (ii > 1) then
            ja(ptr) = kk - 1
            a(ptr)  = aL
            ptr = ptr + 1
        end if

        ! Right
        if (ii < nx) then
            ja(ptr) = kk + 1
            a(ptr)  = aR
            ptr = ptr + 1
        end if

        ! Bottom
        if (jj > 1) then
            ja(ptr) = kk - nx
            a(ptr)  = aB
            ptr = ptr + 1
        end if

        ! Top
        if (jj < ny) then
            ja(ptr) = kk + nx
            a(ptr)  = aT
            ptr = ptr + 1
        end if

        ia(kk+1) = ptr

      end do
    end do

    ia(N+1) = ptr

    ! Trim storage
    a  = a(1:ptr-1)
    ja = ja(1:ptr-1)

  end subroutine assemble_2D_diffusion_matrix

  subroutine assemble_3D_diffusion_matrix(nx, ny, nz, x_domain, y_domain, z_domain, D, Sigma_a, gamma, ia, ja, a, dx, dy, dz)
    use m_matrix_check
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: x_domain, y_domain, z_domain, D, Sigma_a, gamma
    real(8), intent(out) :: dx, dy, dz
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:)

    integer :: N, ii, jj, kk, row, ptr, est_nnz
    integer :: idxL, idxR, idxB, idxT, idxBk, idxFr
    real(8) :: aL, aR, aB, aT, aBack, aFront, aC

    logical :: is_spd, is_symmetric, is_pd

    ! grid spacing
    dx = x_domain / real(nx)
    dy = y_domain / real(ny)
    dz = z_domain / real(nz)

    ! number of unknowns
    N = nx * ny * nz

    ! each row has ≤ 7 nonzeros
    est_nnz = 7*N
    allocate(ia(N+1), ja(est_nnz), a(est_nnz))

    ptr = 1
    ia(1) = 1

    ! loop in k-jj-ii order (k slowest, ii fastest)
    do kk = 1, nz
       do jj = 1, ny
          do ii = 1, nx

             ! flatten to row index
             row = ii + (jj-1)*nx + (kk-1)*nx*ny

             ! stencil coefficients (base values)
             aL    = -D * dy * dz / dx
             aR    = -D * dy * dz / dx
             aB    = -D * dx * dz / dy
             aT    = -D * dx * dz / dy
             aBack = -D * dx * dy / dz
             aFront= -D * dx * dy / dz

             ! diagonal
             aC = Sigma_a*dx*dy*dz - (aL + aR + aB + aT + aBack + aFront)

            !-------------------------------------------
            ! Neumann BC 
            !-------------------------------------------

            !  if (ii == 1)  aC = aC + aL
            !  if (ii == nx) aC = aC + aR
            !  if (jj == 1)  aC = aC + aB
            !  if (jj == ny) aC = aC + aT
            !  if (kk == 1)  aC = aC + aBack
            !  if (kk == nz) aC = aC + aFront

             !-------------------------------------------
             ! Dirichlet BC 
             !-------------------------------------------

            !  if (ii == 1)  aC = aC - aL
            !  if (ii == nx) aC = aC - aR
            !  if (jj == 1)  aC = aC - aB
            !  if (jj == ny) aC = aC - aT
            !  if (kk == 1)  aC = aC - aBack
            !  if (kk == nz) aC = aC - aFront

            !-------------------------------------------
            ! Vacuum BC 
            !-------------------------------------------

            !  if (ii == 1)  aC = aC + aL + (2*D*dy*dz)/(4*D+dx)
            !  if (ii == nx) aC = aC + aR + (2*D*dy*dz)/(4*D+dx)
            !  if (jj == 1)  aC = aC + aB + (2*D*dx*dz)/(4*D+dy)
            !  if (jj == ny) aC = aC + aT + (2*D*dx*dz)/(4*D+dy)
            !  if (kk == 1)  aC = aC + aBack + (2*D*dx*dy)/(4*D+dz)
            !  if (kk == nz) aC = aC + aFront + (2*D*dx*dy)/(4*D+dz)

            
            !-------------------------------------------
            ! Albedo BC 
            !-------------------------------------------

             if (ii == 1)  aC = aC + aL + (2*D*gamma*dy*dz)/(2+gamma*dx)
             if (ii == nx) aC = aC + aR + (2*D*gamma*dy*dz)/(2+gamma*dx)
             if (jj == 1)  aC = aC + aB + (2*D*gamma*dx*dz)/(2+gamma*dy)
             if (jj == ny) aC = aC + aT + (2*D*gamma*dx*dz)/(2+gamma*dy)
             if (kk == 1)  aC = aC + aBack + (2*D*gamma*dx*dy)/(2+gamma*dz)
             if (kk == nz) aC = aC + aFront + (2*D*gamma*dx*dy)/(2+gamma*dz)

            !-------------------------------------------
            ! Write center (always valid)
            !-------------------------------------------
            ja(ptr) = row
            a(ptr)  = aC
            ptr = ptr + 1
            
            ! LEFT
            if (ii > 1 .and. nx > 1) then
              idxL = row - 1
              ja(ptr) = idxL
              a(ptr)  = aL
              ptr = ptr + 1
            end if

            ! RIGHT
            if (ii < nx .and. nx > 1) then
              idxR = row + 1
              ja(ptr) = idxR
              a(ptr)  = aR
              ptr = ptr + 1
            end if

            ! BOTTOM
            if (jj > 1 .and. ny > 1) then
              idxB = row - nx
              ja(ptr) = idxB
              a(ptr)  = aB
              ptr = ptr + 1
            end if

            ! TOP
            if (jj < ny .and. ny > 1) then
              idxT = row + nx
              ja(ptr) = idxT
              a(ptr)  = aT
              ptr = ptr + 1
            end if

            ! BACK (k-1)
            if (kk > 1 .and. nz > 1) then
              idxBk = row - nx*ny
              ja(ptr) = idxBk
              a(ptr)  = aBack
              ptr = ptr + 1
            end if

            ! FRONT (k+1)
            if (kk < nz .and. nz > 1) then
              idxFr = row + nx*ny
              ja(ptr) = idxFr
              a(ptr)  = aFront
              ptr = ptr + 1
            end if

            ! CSR row pointer
            ia(row+1) = ptr

          end do
       end do
    end do

    ia(N+1) = ptr

    ! trim unused (ptr-1 actual nnz)
    ja = ja(1:ptr-1)
    a  = a(1:ptr-1)

    call check_matrix_SPD_CSR(ia, ja, a, is_spd, is_symmetric, is_pd, verbose=.true.)

  end subroutine assemble_3D_diffusion_matrix

  subroutine assemble_rz_diffusion_matrix(nr, nz, r_domain, z_domain, D, Sigma_a, ia, ja, a, dr, dz)
    implicit none
    !---------------------------------------------------------------
    ! Assembles CSR matrix for 2D diffusion operator:
    !   -D*(d2/dx2 + d2/dy2) + Sigma_a
    ! Dirichlet boundary conditions (phi=0) on all edges.
    !
    ! Input:
    !   nx, ny     - grid points (including boundaries)
    !   dx, dy     - grid spacing
    !   D          - diffusion coefficient
    !   Sigma_a    - absorption coefficient
    !
    ! Output:
    !   ia(:)      - kk pointer (CSR)
    !   ja(:)      - column indices (CSR)
    !   a(:)       - nonzero values (CSR)
    !
    ! Notes:
    !   Lexicographic ordering: k = ii + (jj-1)*nx
    !---------------------------------------------------------------
    integer, intent(in) :: nr, nz
    real(8), intent(in) :: r_domain, z_domain, D, Sigma_a
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:)
    real(8), intent(out) ::  dr, dz

    integer :: N, ii, jj, kk, ptr, est_nnz
    real(8) :: aL, aB, aR, aT, aC
    real(8), allocatable :: ri(:)

    dr = r_domain/real(nr) ; dz = z_domain/real(nz) !Linearly spaced dx, dy
    
    N = nr * nz
    est_nnz = 5 * N
    allocate(ia(N+1), ja(est_nnz), a(est_nnz), ri(nr))

    ptr = 1
    ia(1) = 1

    do ii = 1,nr
      ri(ii) = (ii-0.5d0)*dr
    end do

    do jj = 1, nz
      do ii = 1, nr

        kk = ii + (jj-1)*nr

        aL = -D*(ri(ii)-dr/2)/(dr**2*ri(ii)) !-D / dr**2 + D/(2*(ri(ii)-0.5*ii*dr)*dr)
        aR = -D*(ri(ii)+dr/2)/(dr**2*ri(ii)) !-D / dr**2 - D/(2*(ri(ii)+0.5*ii*dr)*dr) 
        aB = -D / (dz**2)
        aT = -D / (dz**2)
        
        ! Multiply absorption by cell volume in r-z (finite-volume form):
        aC = Sigma_a - (aL + aR + aB + aT)
        
        !------------------------------------------------------------
        ! Reflective at r = 0  (Neumann: ∂φ/∂r = 0)
        !------------------------------------------------------------
        if (ii == 1) then
          ! Reflective at r=0 (Neumann ∂φ/∂r = 0): remove the missing left-face
          ! contribution by adding its value to the diagonal (same pattern as 2D/rth)
          aC = aC - aL
          aL = 0.0d0
        end if

        !------------------------------------------------------------
        ! Dirichlet at r = r_max
        !------------------------------------------------------------
        if (ii == nr) then
            aC = aC - aR
            aR = 0.0d0
        end if

        !------------------------------------------------------------
        ! Dirichlet at z = 0
        !------------------------------------------------------------
        if (jj == 1) then
            aC = aC - aB
            aB = 0.0d0
        end if

        !------------------------------------------------------------
        ! Dirichlet at z = z_max
        !------------------------------------------------------------
        if (jj == nz) then
          aC = aC - aT
          aT = 0.0d0
        end if

        ! Center
        ja(ptr) = kk
        a(ptr) = aC
        ptr = ptr + 1

        ! Left
        if (ii > 1) then
          ja(ptr) = kk-1
          a(ptr)  = aL
          ptr = ptr + 1
        endif
        ! Right
        if (ii < nr) then
          ja(ptr) = kk+1
          a(ptr)  = aR
          ptr = ptr + 1
        endif
        ! Down
        if (jj > 1) then
          ja(ptr) = kk-nr
          a(ptr)  = aB
          ptr = ptr + 1
        endif
        ! Up
        if (jj < nz) then
          ja(ptr) = kk+nr
          a(ptr)  = aT
          ptr = ptr + 1
        endif
        ia(kk+1) = ptr
      end do
    end do

    ! Finalize CSR structure
    ia(N+1) = ptr
    if (ptr-1 < size(a)) then
      a = a(1:ptr-1)
      ja = ja(1:ptr-1)
    end if

  end subroutine assemble_rz_diffusion_matrix

  subroutine assemble_rth_diffusion_matrix(nr, nth, ri, dr, dth, D, Sigma_a, ia, ja, a)
    implicit none

    integer, intent(in) :: nr, nth
    real(8), intent(in) :: ri(:), dr, dth, D, Sigma_a
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:)

    integer :: ii, jj, row, nnz, est_nnz, col
    real(8) :: r, aN, aS, aE, aW, a0

    est_nnz = 5 * nr * nth          ! 5-point stencil
    allocate(ia(nr*nth+1), ja(est_nnz), a(est_nnz))

    nnz = 1
    ia(1) = 1

    do ii = 1, nr
      do jj = 1, nth

        row = (ii-1)*nth + jj
        r = ri(ii)

        ! angular diffusion
        aE = -D / (r*r * dth*dth)
        aW =  aE

        ! radial diffusion (symmetric discretization)
        aN = -D * (r+0.5d0*dr) / (r * dr*dr)
        aS = -D * (r-0.5d0*dr) / (r * dr*dr)

        !--- MODIFY FOR BOUNDARIES ---!

        ! Reflective at r=0
        if (ii == 1) then
            ! Neumann BC: dφ/dr=0 → S-coefficient mirrors N
            aS = aN
        end if

        ! Dirichlet at r=R
        if (ii == nr) then
            ! Replace row with φ = 0
            ja(nnz) = row
            a(nnz) = 1.0d0
            nnz = nnz + 1
            ia(row+1) = nnz
            cycle
        end if

        ! center coefficient
        a0 = -(aN + aS + aE + aW) + Sigma_a * r

        !--------- SOUTH ----------
        if (ii > 1) then
            col = (ii-2)*nth + jj
            ja(nnz) = col
            a(nnz)  = aS
            nnz = nnz + 1
        end if

        !---------- NORTH ---------
        col = ii*nth + jj
        ja(nnz) = col
        a(nnz)  = aN
        nnz = nnz + 1

        !---------- WEST ----------
        col = (ii-1)*nth + merge(jj-1, nth, jj>1)
        ja(nnz) = col
        a(nnz)  = aW
        nnz = nnz + 1

        !---------- EAST ----------
        col = (ii-1)*nth + merge(jj+1, 1, jj<nth)
        ja(nnz) = col
        a(nnz)  = aE
        nnz = nnz + 1

        !---------- CENTER ----------
        ja(nnz) = row
        a(nnz)  = a0
        nnz = nnz + 1

        ia(row+1) = nnz

      end do
    end do

    ! trim storage
    a  = a(1:nnz-1)
    ja = ja(1:nnz-1)

  end subroutine

subroutine assemble_3D_diffusion_matrix_g(nx, ny, nz, N, G, gg, x_domain, y_domain, z_domain, D, Sigma_r, gamma, A, dx, dy, dz)
    use CSR_types, only: CSRMatrix
    implicit none

    ! INPUTS
    integer, intent(in) :: nx, ny, nz, N, G, gg
    real(8), intent(in) :: x_domain, y_domain, z_domain
    real(8), intent(in) :: D, Sigma_r, gamma

    ! OUTPUTS
    real(8), intent(out) :: dx, dy, dz
    type(CSRMatrix), intent(out) :: A   ! output sparse matrix

    ! locals
    integer :: ii, jj, kk, row, ptr, est_nnz
    integer :: idxL, idxR, idxB, idxT, idxBk, idxFr
    real(8) :: aL, aR, aB, aT, aBack, aFront, aC

    ! grid spacing
    dx = x_domain / real(nx)
    dy = y_domain / real(ny)
    dz = z_domain / real(nz)

    ! estimate nnz (7-point stencil maximum)
    est_nnz = 7 * N

    ! allocate CSR arrays
    allocate(A%ia(N+1))
    allocate(A%ja(est_nnz))
    allocate(A%aa(est_nnz))

    ! initialize
    A%ia = 0
    A%ja = 0
    A%aa = 0

    ptr = 1
    A%ia(1) = 1  ! 1-based CSR convention: IA(1) = index of first element = 1

    ! loop rows in k-j-i order (consistent with your flattening)
    do kk = 1, nz
        do jj = 1, ny
            do ii = 1, nx

                ! flatten to row index (1-based)
                row = ii + (jj-1)*nx + (kk-1)*nx*ny

                ! stencil coefficients (off-diagonals)
                aL = -D * ( dy * dz ) / dx
                aR = aL
                aB = -D * ( dx * dz ) / dy
                aT = aB
                aBack  = -D * ( dx * dy ) / dz
                aFront = aBack

                ! diagonal contribution (positive) - note sign convention: off-diagonals are negative
                aC = Sigma_r * dx * dy * dz - (aL + aR + aB + aT + aBack + aFront)

                ! apply Albedo BC by shifting off-diagonal contributions into diagonal (keeps same sign convention)
                if (ii == 1)  aC = aC + aL + ( 2.0d0*D*gamma*dy*dz )/( 2.0d0 + gamma*dx )
                if (ii == nx) aC = aC + aR + ( 2.0d0*D*gamma*dy*dz )/( 2.0d0 + gamma*dx )
                if (jj == 1)  aC = aC + aB + ( 2.0d0*D*gamma*dx*dz )/( 2.0d0 + gamma*dy )
                if (jj == ny) aC = aC + aT + ( 2.0d0*D*gamma*dx*dz )/( 2.0d0 + gamma*dy )
                if (kk == 1)  aC = aC + aBack + ( 2.0d0*D*gamma*dx*dy )/( 2.0d0 + gamma*dz )
                if (kk == nz) aC = aC + aFront + ( 2.0d0*D*gamma*dx*dy )/( 2.0d0 + gamma*dz )

                ! --- center (column = row) ---
                if (ptr > est_nnz) then
                    print *, "ERROR: ptr exceeds est_nnz in assemble_3D_diffusion_matrix_g, increase est_nnz"
                    stop
                end if
                A%ja(ptr) = row       ! column index of diagonal
                A%aa(ptr) = aC
                ptr = ptr + 1

                ! LEFT neighbour (col = row-1)
                if (ii > 1) then
                    idxL = row - 1
                    A%ja(ptr) = idxL
                    A%aa(ptr) = aL
                    ptr = ptr + 1
                end if

                ! RIGHT neighbour
                if (ii < nx) then
                    idxR = row + 1
                    A%ja(ptr) = idxR
                    A%aa(ptr) = aR
                    ptr = ptr + 1
                end if

                ! BOTTOM neighbour (j-1)
                if (jj > 1) then
                    idxB = row - nx
                    A%ja(ptr) = idxB
                    A%aa(ptr) = aB
                    ptr = ptr + 1
                end if

                ! TOP neighbour (j+1)
                if (jj < ny) then
                    idxT = row + nx
                    A%ja(ptr) = idxT
                    A%aa(ptr) = aT
                    ptr = ptr + 1
                end if

                ! BACK neighbour (k-1)
                if (kk > 1) then
                    idxBk = row - nx*ny
                    A%ja(ptr) = idxBk
                    A%aa(ptr) = aBack
                    ptr = ptr + 1
                end if

                ! FRONT neighbour (k+1)
                if (kk < nz) then
                    idxFr = row + nx*ny
                    A%ja(ptr) = idxFr
                    A%aa(ptr) = aFront
                    ptr = ptr + 1
                end if

                ! set IA for next row (first index after this row's entries)
                A%ia(row+1) = ptr

            end do
        end do
    end do

    ! final NNZ index pointer: ptr is (nnz + 1)
    A%ia(N+1) = ptr

    ! trim JA and AA to actual nnz
    if (ptr - 1 < size(A%ja)) A%ja = A%ja(1:ptr-1)
    if (ptr - 1 < size(A%aa)) A%aa = A%aa(1:ptr-1)

    return
end subroutine assemble_3D_diffusion_matrix_g

end module m_diffusion_matrix

