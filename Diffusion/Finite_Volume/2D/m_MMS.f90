
!------------------------------------------------------------------------!
! Purpose:                                                              -!
!                                                                       -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 19/11/25     T. Charlton    Implemented MMS diffusion                 -! 
!------------------------------------------------------------------------! 

module m_MMS
  use m_constants
  implicit none
  contains
  subroutine MMS_2D_xy_diffusion(dx, dy, nx, ny, x_domain, y_domain, phi)
      integer, intent(in) :: nx, ny
      real(8), intent(in) :: phi(:), dx, dy, x_domain, y_domain

      integer :: ii, j_mid, N
      real(8) :: y_mid, phi_exact, L2_error

      N = nx*ny
      j_mid = (ny + 1) / 2
      y_mid = (j_mid - 0.5d0) * dy

      do ii = 1, nx
          phi_exact = sin(pi*(ii - 0.5d0) * dx/x_domain) * sin(pi*y_mid/y_domain)
          L2_error = L2_error + (phi((j_mid-1)*nx + ii) - phi_exact)**2
      end do
      L2_error = sqrt(L2_error / nx)
      
      open(unit=991, file="flux.dat", status="replace", action="write")

      do ii = 1, nx
          phi_exact = sin(pi*(ii - 0.5d0) * dx/x_domain) * sin(pi*y_mid/y_domain)
          write(991,'(4(1X,E15.8))') ii*dx, phi((j_mid-1)*nx + ii), phi_exact,  L2_error
      end do

      close(991)
  end subroutine MMS_2D_xy_diffusion

end module