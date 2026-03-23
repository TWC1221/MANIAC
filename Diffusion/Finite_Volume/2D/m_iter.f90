
!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  Conduct 2D NDE eigenvalue problems in xy, rz and rth geometries      -!  
!  Input data defining problem geometry, discretization and meshing     -!
!  Input data defining neutron/material interaction                     -!
!  Call functions to assemble & solve coefficient matrix & output VTK   -!
!    - assemble_XX_diffusion_matrix                                     -!
!    - PCG_algorithm                                                    -!
!    - outpVTK                                                          -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 12/11/25     T. Charlton    Implemented 2D matrix                     -!
! 15/11/25     T. Charlton    Implemented rth matrix                    -!
! 17/11/25     T. Charlton    Implemented rz matrix                     -! 
!------------------------------------------------------------------------! 

module m_iter
  use m_diffusion_matrix
  use m_PCG_solver
  use m_outpVTK
  use m_MMS
  use m_constants
  implicit none
  contains

  subroutine iter(PCG_mode, max_iter)
    implicit none

    integer, intent(in) :: PCG_mode, max_iter
    
    integer, allocatable :: ia(:), ja(:)
    real(8), allocatable :: aa(:), phi(:), phi1(:), lambda(:), src(:)
    real(8) :: norm, dx, dy
    integer :: N, ii = 1 !jj = 1

    real(8) :: D = 1/0.6
    real(8) :: Sigma_a = 0.2
    real(8) :: nuSigma_f = 1.0

    integer :: nx = 11, ny = 12
    real(8) :: x_domain = 20.0
    real(8) :: y_domain = 20.0

    call system("pkill gnuplot")

    N = nx*ny
    allocate(phi(N), phi1(N), lambda(10000), src(N))
    
    call assemble_2D_diffusion_matrix(nx, ny, x_domain, y_domain, D, Sigma_a, ia, ja, aa, dx, dy)

    phi = 1
    lambda(1) = 1

    ! do ii = 1, nx
    !   do jj = 1, ny
    !     src((jj-1)*nx+ii) = sin(pi*dx*(ii-0.5d0)/x_domain)*sin(pi*dy*(jj-0.5d0)/y_domain)*(Sigma_a+D*((pi/x_domain)**2+(pi/y_domain)**2))*dx*dy
    !   end do
    ! end do

    do
      ii = ii + 1

      src = (nuSigma_f * phi * dx * dy) / lambda(ii-1)

      phi1 = phi(1:N)
      call PCG_algorithm(aa, ja, ia, phi1, src, PCG_mode, max_iter)

      ! call MMS_2D_xy_diffusion(dx, dy, nx, ny, x_domain, y_domain, phi1)
      ! call outpVTK(phi1, nx, ny, dx, dy)
      ! stop

      lambda(ii) = lambda(ii-1) * sum(nuSigma_f*phi1*dx*dy)/sum(nuSigma_f*phi*dx*dy)

      norm = (sum(phi1*nuSigma_f*dx*dy))
      phi = phi1/norm

      print*,'Iteration =', ii, 'Keff =', lambda(ii), 'norm_check =', sum(phi*nuSigma_f*dx*dy)/(lambda(ii))

      if (abs((lambda(ii) - lambda(ii-1)) / lambda(ii-1)) < 1.0d-13) exit
      if (maxval(abs(phi1 - phi)) < 1.0d-12) exit
    end do

    call outpVTK(phi1, nx, ny, dx, dy)

  end subroutine iter

  subroutine iter_xyz(PCG_mode, max_iter)
    implicit none

    integer, intent(in) :: PCG_mode, max_iter
    
    integer, allocatable :: ia(:), ja(:)
    real(8), allocatable :: aa(:), phi(:), phi1(:), lambda(:), src(:), src_ext(:)
    real(8) :: alpha, gamma, norm, dx, dy, dz
    integer :: N, ii, jj, kk, row 
    ! NEEDED FOR SOURCE SURFACE

    real(8) :: D = 1.0/0.6
    real(8) :: Sigma_a = 0.2
    real(8) :: nuSigma_f = 1.0 !1.0

    integer :: nx = 90
    integer :: ny = 90
    integer :: nz = 90
    real(8) :: x_domain = 20.0
    real(8) :: y_domain = 20.0
    real(8) :: z_domain = 20.0

    alpha = 0.0
    gamma = (1-alpha)/(2*D*(1+alpha))

    call system("pkill gnuplot")

    N = nx*ny*nz
    allocate(phi(N), phi1(N), lambda(max_iter), src(N),src_ext(N))
    
    call assemble_3D_diffusion_matrix(nx, ny, nz, x_domain, y_domain, z_domain, D, Sigma_a, gamma, ia, ja, aa, dx, dy, dz)

    !Source RHS term
    ! do kk = 1, nz
    !   do jj = 1, ny
    !     do ii = 1, 1
    !       row = ii + (jj-1)*nx + (kk-1)*nx*ny
    !       src_ext(row) = D*gamma*10/(2+gamma*dx)
    !     end do
    !   end do 
    ! end do

    ii = 1
    phi = 1
    lambda(1) = 1

    do
      ii = ii + 1

      src = ((nuSigma_f * phi) * dx*dy*dz) / lambda(ii-1) ! + src_ext

      phi1 = phi
      call PCG_algorithm(aa, ja, ia, phi1, src, PCG_mode, max_iter)

      !exit !Source Term run path (comment out lambda below)
      lambda(ii) = lambda(ii-1) * ( sum(phi1 * nuSigma_f * dx*dy*dz) / sum(phi * nuSigma_f * dx*dy*dz) )

      norm = sum(phi1 * nuSigma_f * dx*dy*dz) / lambda(ii)
      phi  = phi1 / norm

      print*,'Iteration =', ii, 'Keff =', lambda(ii), 'norm_check =', sum(phi*nuSigma_f * dx*dy*dz)/(lambda(ii))

      if (abs((lambda(ii) - lambda(ii-1)) / lambda(ii-1)) < 1.0d-10) exit
      !if (maxval(abs(phi1 - phi)) < 1.0d-13) then; print*,('phi'); exit; end if
    end do

    call outpVTK_xyz(phi, nx, ny, nz, dx, dy, dz) !phi 1 for source term

  end subroutine iter_xyz

  subroutine iter_rz(PCG_mode, max_iter)
    implicit none
    integer, intent(in) :: PCG_mode, max_iter
    integer, allocatable :: ia(:), ja(:)
    real(8), allocatable :: aa(:), phi(:), phi1(:), phi_ptr(:), lambda(:), src(:), r_node(:)
    real(8) :: dr, dz
    integer :: nr = 100, nz = 100, N, ii
    !real(8) :: r_domain = 2.0, z_domain = 2.0, D = 0.04, Sigma_a = 0.02, nuSigma_f = 3.0, norm
    real(8) :: r_domain = 20.0, z_domain = 20.0, D = 1/0.6, Sigma_a = 0.2, nuSigma_f = 1.0, norm

    call system("pkill gnuplot")

    N = nr*nz
    allocate(phi(N), phi1(N), phi_ptr(N), lambda(100000), src(N), r_node(N))

    call assemble_rz_diffusion_matrix(nr, nz, r_domain, z_domain, D, Sigma_a, ia, ja, aa, dr, dz)
    
    phi(1:N) = 1; lambda(1) = 1; ii = 1

    do
      ii = ii + 1
      ! Scale source by cell area in rz geometry
      src = (nuSigma_f * phi) / lambda(ii-1)

      phi_ptr = phi(1:N)
      call PCG_algorithm(aa, ja, ia, phi_ptr, src, PCG_mode, max_iter)
      phi1(1:N) = phi_ptr

      lambda(ii) = lambda(ii-1) * sum(nuSigma_f*phi1*dr*dz)/sum(nuSigma_f*phi*dr*dz)

      norm = sum(phi1*nuSigma_f*dr*dz)
      phi = phi1/norm

      !print*,phi1(10),phi(10),lambda(ii), lambda(ii-1)
      print*,"Iteration =", ii, "Keff =", lambda(ii), "norm_check =",sum(phi*nuSigma_f*dr*dz)/(lambda(ii))

      if (abs((lambda(ii) - lambda(ii-1)) / lambda(ii-1)) < 1.0d-8) exit
    end do

    !print'(20F6.2)',phi1
    call outpVTK_rz(phi, nr, nz, dr, dz)

  end subroutine iter_rz

! subroutine iter_rth(PCG_mode, max_pcg)
!     implicit none
!     !-----------------------------
!     ! Inputs
!     !-----------------------------
!     integer, intent(in) :: PCG_mode, max_pcg

!     !-----------------------------
!     ! Local variables
!     !-----------------------------
!     integer :: nr, nth, N, i, j, iter, max_outer
!     real(8) :: R_domain, D, Sigma_a, nuSigma_f
!     real(8) :: dr, dth
!     real(8) :: k_old, k_new, rel_change
!     real(8) :: vol_weight

!     !-----------------------------
!     ! Arrays
!     !-----------------------------
!     integer, allocatable :: ia(:), ja(:)
!     real(8), allocatable :: aa(:)
!     real(8), allocatable :: phi(:), phi_new(:), src(:)
!     real(8), allocatable :: r(:), r_node(:)

!     !--------------------------------------
!     ! Geometry and physics
!     !--------------------------------------
!     nr        = 100
!     nth       = 200
!     N         = nr*nth
!     R_domain  = 10d0

!     D         = 1d0/0.6d0
!     Sigma_a   = 0.2d0
!     nuSigma_f = 0.01d0

!     max_outer = 200

!     !--------------------------------------
!     ! Allocate arrays
!     !--------------------------------------
!     allocate(phi(N), phi_new(N), src(N))
!     allocate(r(nr), r_node(N))

!     !--------------------------------------
!     ! Radial coordinates (cell centers)
!     !--------------------------------------
!     dr  = R_domain/real(nr)
!     dth = 2.0d0*pi/real(nth)

!     do i = 1, nr
!         r(i) = (i-0.5d0)*dr
!         do j = 1, nth
!             r_node((i-1)*nth + j) = r(i)
!         end do
!     end do

!     !--------------------------------------
!     ! Assemble CSR diffusion matrix
!     !--------------------------------------
!     call assemble_rth_diffusion_matrix(nr, nth, r, dr, dth, D, Sigma_a, ia, ja, aa)
!     print *, "Matrix Assembly Successful"

!     !--------------------------------------
!     ! Initial guess
!     !--------------------------------------
!     phi = 1.0d0
!     k_old = 1.0d0

!     !--------------------------------------
!     ! Power iteration (flux solve + k update)
!     !--------------------------------------
!     do iter = 1, max_outer

!         !--- Build fixed source: fission / keff_old
!         src = (nuSigma_f * phi * r_node) / k_old

!         ! PCG solve:  A * phi_new = src
!         phi_new = phi
!         call PCG_algorithm(aa, ja, ia, phi_new, src, PCG_mode, max_pcg)

!         !--------------------------------------
!         ! Update keff using the Rayleigh quotient:
!         !
!         !  k_new = k_old * ( ∫ νΣf φ_new dV ) / ( ∫ νΣf φ dV )
!         !--------------------------------------
!         k_new = k_old * sum(nuSigma_f * phi_new * r_node) / &
!                          sum(nuSigma_f * phi     * r_node)

!         !--------------------------------------
!         ! Normalize flux (L2 weighted by r dV)
!         !--------------------------------------
!         vol_weight = sqrt(sum((phi_new**2) * r_node * dr * dth))
!         phi_new = phi_new / vol_weight

!         !--------------------------------------
!         ! Print diagnostics
!         !--------------------------------------
!         print '(A,I4,A,ES20.10,A,ES20.10)', &
!             "Iter=", iter, "   keff=", k_new, "   norm_check=", &
!             sum(nuSigma_f * phi_new * r_node * dr * dth) / k_new

!         !--------------------------------------
!         ! Convergence test (relative keff change)
!         !--------------------------------------
!         rel_change = abs((k_new - k_old)/k_old)

!         if (rel_change < 1d-10) exit

!         ! Prepare for next iteration
!         phi = phi_new
!         k_old = k_new

!     end do

!     !--------------------------------------
!     ! Output
!     !--------------------------------------
!     call outpVTK_rth(phi_new, nr, nth, dr, dth)

! end subroutine iter_rth


end module