module m_MMS_solver
!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Contains the subroutines required to solve the tridiagonal matrix    -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 08/10/2025    T. Charlton      Original code                          -!
!------------------------------------------------------------------------! 

implicit none
contains
    subroutine build_matrix_A_MMS(a, b, c, N, phi, L2Err)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N
        real(8), intent(inout) :: a(N), b(N), c(N)
        real(8), intent(out) :: phi(N), L2Err

        real(8), parameter :: L = 1.0
        real(8), parameter :: Dif = 1.0/3.0
        real(8), parameter :: Sigma_a = 0.01
        real(8), parameter :: S0 = 1.0
        real(8) :: dx, phi_exact(N)
        real(8), dimension(N) :: d
        integer :: i

        dx = L/(N-1)

        do i = 1, N
            if (i == 1) then
                a(i) = 0.0
                b(i) = VLARGE_NUMBER !Sigma_a + 2*Dif/dx**2 + (1/dx)*(1-alpha)/(1+alpha)
                c(i) = -2*Dif/dx**2
                d(i) = (pi**2*Dif + Sigma_a) * sin(pi*(i-1)*dx)
            elseif (i == N) then
                a(i) = -2*Dif/dx**2
                b(i) = VLARGE_NUMBER !Sigma_a + 2.0*Dif/dx**2 + (1/dx)*(1-alpha)/(1+alpha)
                c(i) = 0.0
                d(i) = (pi**2*Dif + Sigma_a) * sin(pi*(i-1)*dx)
            else
                a(i) = -Dif/dx**2
                b(i) = Sigma_a + 2.0*Dif/dx**2
                c(i) = -Dif/dx**2
                d(i) = (pi**2*Dif + Sigma_a) * sin(pi*(i-1)*dx)
            end if
        end do

        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

        ! Write all columns at once after the loop
        open(unit=991, file="flux.dat", status="replace", action="write")

        do i = 1, N
            write(991,'(F10.5,2(1X,E15.8))') (i-1)*dx, phi(i), sin(pi*(i-1)*dx)
            phi_exact(i) = sin(pi*(i-1)*dx)
        end do
        close(991)

        L2Err = sqrt(sum((phi - phi_exact)**2)/N)

         !call inverse_power_iteration(a, b, c, N, max_iter, tol, lambda, phi)

    end subroutine build_matrix_A_MMS

end module m_MMS_solver
