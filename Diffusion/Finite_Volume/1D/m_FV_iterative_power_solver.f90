module m_iterative_power_solver
!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Contains the subroutines required to solve the tridiagonal matrix    -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 08/10/2025    T. Charlton      Original code                          -!
! 10/11/2025    T. Charlton      Finite Volume Implentation             -!
!------------------------------------------------------------------------! 

implicit none
contains
    subroutine build_matrix_A_iterative_power(a, b, c, N, alpha, phi, dx, Dif, Sigma_a, S0)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: alpha
        real(8), intent(inout) :: a(N), b(N), c(N), phi(N), dx(N), Dif(N), Sigma_a(N), S0(N)

        real(8), dimension(N) :: d, x, Qsurf
        integer :: ii

        Qsurf = 10
        do ii = 1, N
            x(ii) =  sum(dx(1:ii))
        end do

        do ii = 1, N
            if (ii == 1) then
                a(ii) = 0.0d0
                ! Albedo
                !b(ii) =  2*(Dif(ii+1)*Dif(ii))/(dx(ii) * (Dif(ii+1)*dx(ii) + Dif(ii)*dx(ii+1))) + 2*(Dif(ii)*(1-alpha))/(dx(ii)*(4*Dif(ii)*(1+alpha)+dx(ii)*(1-alpha))) + Sigma_a(ii)
                ! Zero flux
                !b(ii) =  2*(Dif(ii+1)*Dif(ii))/(dx(ii) * (Dif(ii+1)*dx(ii) + Dif(ii)*dx(ii+1))) + 2*Dif(ii)/dx(ii)**2 + Sigma_a(ii)
                ! Surface source
                b(ii) =  2*(Dif(ii+1)*Dif(ii))/(dx(ii) * (Dif(ii+1)*dx(ii) + Dif(ii)*dx(ii+1))) + 2*Dif(ii)/(dx(ii)*(4*Dif(ii)+dx(ii))) + Sigma_a(ii)
                c(ii) = -2.0d0*(Dif(ii+1)*Dif(ii))/(dx(ii) * (Dif(ii+1)*dx(ii) + Dif(ii)*dx(ii+1)))
                ! Source
                !d(ii) = S0(ii)
                ! Surface Source
                d(ii) = S0(ii)+4*Qsurf(ii)*Dif(ii)/(dx(ii)*(4*Dif(ii)+dx(ii)))
            elseif (ii == N) then
                a(ii) = -2*(Dif(ii-1)*Dif(ii))/(dx(ii-1) * (Dif(ii-1)*dx(ii) + Dif(ii)*dx(ii-1)))
                ! Albedo
                !b(ii) =  2*(Dif(ii-1)*Dif(ii))/(dx(ii-1) * (Dif(ii-1)*dx(ii) + Dif(ii)*dx(ii-1))) + 2*(Dif(ii)*(1-alpha))/(dx(ii)*(4*Dif(ii)*(1+alpha)+dx(ii)*(1-alpha))) + Sigma_a(ii)
                ! Zero flux
                !b(ii) =  2*(Dif(ii-1)*Dif(ii))/(dx(ii-1) * (Dif(ii-1)*dx(ii) + Dif(ii)*dx(ii-1))) + 2*Dif(ii)/dx(ii)**2 + Sigma_a(ii)
                ! Surface source
                b(ii) =  2*(Dif(ii-1)*Dif(ii))/(dx(ii-1) * (Dif(ii-1)*dx(ii) + Dif(ii)*dx(ii-1))) + 2*Dif(ii)/(dx(ii)*(4*Dif(ii)+dx(ii))) + Sigma_a(ii)
                c(ii) = 0.0d0
                ! Source
                !d(ii) = S0(ii)
                ! Surface Source
                d(ii) = S0(ii)+4*Qsurf(ii)*Dif(ii)/(dx(ii)*(4*Dif(ii)+dx(ii)))
            else
                a(ii) = -2*(Dif(ii-1)*Dif(ii))/(dx(ii-1) * (Dif(ii-1)*dx(ii) + Dif(ii)*dx(ii-1)))
                b(ii) = 2*(Dif(ii-1)*Dif(ii))/(dx(ii)*(Dif(ii-1)*dx(ii) + Dif(ii)*dx(ii-1))) + 2*(Dif(ii+1)*Dif(ii))/(dx(ii)*(Dif(ii+1)*dx(ii) + Dif(ii)*dx(ii+1))) +Sigma_a(ii) 
                c(ii) = -2*(Dif(ii+1)*Dif(ii))/(dx(ii) * (Dif(ii+1)*dx(ii) + Dif(ii)*dx(ii+1)))
                d(ii) = S0(ii)
            end if
        end do
        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

        ! Write all columns at once after the loop
        open(unit=991, file="flux.dat", status="replace", action="write")
        do ii = 1, N
            write(991,'(F10.5,2(1X,E15.8))') x(ii), phi(ii)
        end do
        close(991)
    end subroutine build_matrix_A_iterative_power

end module m_iterative_power_solver
