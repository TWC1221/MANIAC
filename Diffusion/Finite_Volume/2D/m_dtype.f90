!--------------------------------------------------------------------------------------------!
!  Purpose:                                                                                 -!
!  Contains inner/outer iteration on multigroup diffusion eigenvalue problem                -!
!           no-upscatter fixed fission source solver                                        -! 
!--------------------------------------------------------------------------------------------!
!  Record of revisions:                                                                     -!
!   Date       Programmer     Description of change                                         -!
!   ====       ==========     =====================                                         -!
! 27/11/2025    T. W Charlton     Original code                                             -!
!--------------------------------------------------------------------------------------------! 
module multigroup_multidimensional
    USE m_PCG_solver
    use m_diffusion_matrix
    implicit none
    contains

subroutine multigroup_diffusion_iter_3d(PCG_mode, max_iter)
    use CSR_types, only: CSRMatrix
    implicit none

    integer, intent(in) :: PCG_mode, max_iter

    ! local problem size
    integer :: G, nx, ny, nz, N
    integer :: ii, jj, gg, i, j
    real(8) :: X_domain, Y_domain, Z_domain
    real(8) :: dx, dy, dz, alpha

    ! CSR matrices for each group
    type(CSRMatrix), allocatable :: A_all(:)

    ! fields and cross-sections
    real(8), allocatable :: phi(:,:), phi_prime(:,:), phi_ptr(:)
    real(8), allocatable :: D(:), Sigma_t(:), Sigma_a(:), Sigma_r(:), nu_Sigma_f(:)
    real(8), allocatable :: Sigma_s(:,:), Sigma_s_upscatter(:,:), Sigma_s_L(:,:), Sigma_s_U(:,:)
    real(8), allocatable :: Sigma_s_L_sum(:,:), Sigma_s_U_sum(:,:)
    real(8), allocatable :: S_f(:), S_f_iter(:), K_eff(:), x(:), chi(:)

    !------------------------------
    ! Problem size and groups
    !------------------------------
    G = 3
    nx = 40
    ny = 40
    nz = 40
    N = nx * ny * nz

    X_domain = 10.0d0
    Y_domain = 10.0d0
    Z_domain = 10.0d0
    alpha = 0.0d0

    dx = X_domain / nx
    dy = Y_domain / ny
    dz = Z_domain / nz

    !------------------------------
    ! allocate arrays
    !------------------------------
    allocate(A_all(G))
    allocate(phi(G, N))
    allocate(phi_prime(G, N))
    allocate(phi_ptr(N))
    allocate(D(G), Sigma_t(G), Sigma_a(G), Sigma_r(G), nu_Sigma_f(G))
    allocate(Sigma_s(G,G), Sigma_s_upscatter(G,G), Sigma_s_L(G,G), Sigma_s_U(G,G))
    allocate(Sigma_s_L_sum(G,N), Sigma_s_U_sum(G,N))
    allocate(S_f(N), S_f_iter(N), K_eff(max_iter+1), x(N), chi(G))

    !------------------------------
    ! Initialize physics
    !------------------------------
    chi = [0.9d0, 0.1d0, 0.0d0]
    Sigma_a = [0.015d0, 0.04d0, 0.12d0]
    nu_Sigma_f = [0.02d0, 0.1d0, 0.35d0]

    ! initial flux guess
    phi = 1.0d-6   ! small positive initial guess to avoid zeros
    phi_prime = 0.0d0

    S_f = 1.0d0    ! initial fission source (will be overwritten)
    K_eff(1) = 1.0d0

    ! Test scattering matrix (example)
    Sigma_s_upscatter = reshape([0.20d0, 0.10d0, 0.17d0, &
                                 0.05d0, 0.25d0, 0.10d0, &
                                 0.10d0, 0.07d0, 0.30d0], shape=[3,3])

    ! split into L (downscatter to lower groups) and U (upscatter to higher groups)
    Sigma_s_L = 0.0d0
    Sigma_s_U = 0.0d0

    do i = 1, G
        do j = 1, G
            if (i > j) then
                Sigma_s_L(i,j) = Sigma_s_upscatter(i,j)  ! i>j: downscatter (to lower group index)
            else if (i < j) then
                Sigma_s_U(i,j) = Sigma_s_upscatter(i,j)  ! i<j: upscatter
            end if
            ! diagonal scattering removed (if present) - keep diagonal zero in L and U
        end do
    end do

    !------------------------------
    ! Assemble diffusion matrices for each group (CSR)
    !------------------------------
    do gg = 1, G
        ! compute total and removal cross sections for diffusion
        Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s_upscatter(gg,1:G))
        Sigma_r(gg) = Sigma_a(gg) + sum(Sigma_s_upscatter(gg,1:G)) - Sigma_s_upscatter(gg,gg)
        D(gg) = 1.0d0 / (3.0d0 * Sigma_t(gg))

        ! call assemble; returns A_all(gg) with AA/JA/IA set
        call assemble_3D_diffusion_matrix_g(nx, ny, nz, N, G, gg, X_domain, Y_domain, Z_domain, D(gg), Sigma_r(gg), alpha, A_all(gg), dx, dy, dz)
    end do

    !------------------------------
    ! Outer power iteration with inner solve for each group (Gauss-Seidel-like)
    !------------------------------
    do ii = 2, max_iter
        S_f_iter = S_f

        ! inner iterative cycle that updates phi for each group until inner convergence
        do
            ! for each energy group update RHS and solve
            do gg = 1, G

                ! Build scattering contributions into RHS:
                ! S_f = chi(gg) / K_eff * previous_fission_source + sum_downscatter(L*phi_prime) + sum_upscatter(U*phi)
                ! compute downscatter (i>j) contributions using phi_prime (previous inner iterate)
                Sigma_s_L_sum(:,:) = 0.0d0
                Sigma_s_U_sum(:,:) = 0.0d0

                do jj = 1, G
                    ! Note: Sigma_s_L(jj,gg) multiplies flux of group jj and scatters into group gg
                    ! So contributions for group gg are Sigma_s_L(gg,jj)*phi_prime(jj,:) if you use (target,source) ordering.
                    ! Adjust indexing to match your stored ordering: here we assume Sigma_s_L(target,source)
                    Sigma_s_L_sum(jj,1:N) = Sigma_s_L(gg,jj) * phi_prime(jj,1:N)
                    Sigma_s_U_sum(jj,1:N) = Sigma_s_U(gg,jj) * phi(jj,1:N)
                end do

                ! sum along source index (jj) to get full scattering contribution into group gg
                ! sum(..., dim=1) yields an array length N
                S_f(1:N) = ((chi(gg) / K_eff(ii-1)) * S_f_iter(1:N) + sum(Sigma_s_L_sum, dim=1) + sum(Sigma_s_U_sum, dim=1))*dx*dy*dz

                ! RHS ready in S_f(1:N); solve A_all(gg) * phi_ptr = S_f using PCG
                phi_ptr = phi(gg, :)
                call PCG_algorithm(A_all(gg)%aa, A_all(gg)%ja, A_all(gg)%ia, phi_ptr, S_f, PCG_mode, max_iter)
                phi(gg, :) = phi_ptr

            end do

            ! ! check inner convergence across groups 2 and 3 (you used these earlier)
            ! if ( maxval( abs( phi(3,:) - phi_prime(3,:) ) ) < 1.0d-10 .and. &
            !      maxval( abs( phi(2,:) - phi_prime(2,:) ) ) < 1.0d-10 ) exit

            ! update phi_prime with current iterate (for downscatter use)
            do gg = 1,G
                if ( maxval( abs( phi(gg,:) - phi_prime(gg,:) ) ) < 1.0d-10) exit
                phi_prime(gg,:) = phi(gg,:)
            end do 

            exit
        end do  ! inner

        ! compute new fission source and keff
        S_f = matmul(transpose(phi), nu_Sigma_f(1:G))
        K_eff(ii) = K_eff(ii-1) * sum(S_f * dx * dy * dz) / sum(S_f_iter * dx * dy * dz)

        print *, "Iteration =", ii-1, "  Keff =", K_eff(ii-1)

        ! check outer convergence (relative keff change and fission source)
        if ( abs( (K_eff(ii) - K_eff(ii-1)) / K_eff(ii-1) ) < 1.0d-10 .and. &
             maxval( abs( S_f - S_f_iter ) ) < 1.0d-10 ) then
            exit
        end if

    end do  ! outer

    do gg = 1,G
        print*,sum(phi(gg,:)*dx*dy*dz)
    end do

    return
end subroutine multigroup_diffusion_iter_3d

end module multigroup_multidimensional

! module vector_module
!     implicit none
!     type :: vector
!         real :: x
!         real :: y
!     end type vector
!     contains
!     type (vector) function vector_add ( v1, v2 )
!         implicit none
!         type (vector), intent(in) :: v1
!         type (vector), intent(in) :: v2

!         vector_add%x = v1%x + v2%x
!         vector_add%y = v1%y + v2%y
!     end function vector_add
! end module vector_module

! PROGRAM test_vectors
!     use vector_module
!     implicit none

!     type(vector) :: v1
!     type(vector) :: v2
!     write (*,*) 'Enter the first vector (x1,y1):'
!     read (*,*) v1%x, v1%y

!     write (*,*) 'Enter the second vector (x2,y2):'
!     read (*,*) v2%x, v2%y

!     write (*,1000) vector_add(v1,v2)
!     1000 FORMAT('The sum of the points (v1, v2) is (',F8.2,','F8.2,')')
! END PROGRAM



