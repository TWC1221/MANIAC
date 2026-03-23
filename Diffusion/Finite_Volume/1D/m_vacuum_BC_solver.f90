module m_Vacuum_BC_solver
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
    subroutine build_matrix_A_vacuum(a, b, c, N, alpha, phi, phi_analytical, L2Err)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: alpha
        real(8), intent(inout) :: a(N), b(N), c(N)
        real(8), intent(out) :: phi(N), L2Err, phi_analytical(N)

        real(8), parameter :: R_domain = 1.0
        real(8), parameter :: Dif = 1.0/3.0
        real(8), parameter :: Sigma_a = 0.01
        real(8), parameter :: S0 = 1.0
        real(8) :: dx, L
        real(8) :: AA
        real(8), dimension(N) :: d
        integer :: ii

        dx = R_domain/(N-1)
        L=sqrt(Dif/Sigma_a)

        do ii = 1, N
            if (ii == 1) then
                a(ii) = 0.0
                b(ii) = Sigma_a + 2*Dif/dx**2 + (1/dx)*(1-alpha)/(1+alpha)
                c(ii) = -2*Dif/dx**2
                d(ii) = S0
            elseif (ii == N) then
                a(ii) = -2*Dif/dx**2
                b(ii) = Sigma_a + 2.0*Dif/dx**2 + (1/dx)*(1-alpha)/(1+alpha)
                c(ii) = 0.0
                d(ii) = S0
            else
                a(ii) = -Dif/dx**2
                b(ii) = Sigma_a + 2.0*Dif/dx**2
                c(ii) = -Dif/dx**2
                d(ii) = S0
            end if
        end do

        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

        ! Write all columns at once after the loop
        open(unit=991, file="flux.dat", status="replace", action="write")

        do ii = 1, N
            AA = (-((S0/(2*Sigma_a))*((L/(2*Dif))*sinh(R_domain/L)+1+cosh(R_domain/L)))/(cosh(R_domain/L)+((L/(4*Dif)+(Dif/L)))*sinh(R_domain/L)))
            phi_analytical(ii) = AA*cosh(dx*(ii-1)/L) + (L/(2*Dif))*(AA + S0/Sigma_a)*sinh(dx*(ii-1)/L) + S0/Sigma_a
            write(991,'(F10.5,2(1X,E15.8))') (ii-1)*dx, phi(ii), phi_analytical(ii)
        end do
        close(991)

        L2Err = sqrt(sum((phi - phi_analytical)**2)/N)

    end subroutine build_matrix_A_vacuum

    subroutine build_matrix_A_heterogeneous_m1(a, b, c, N, N_interface, alpha, phi)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N, N_interface
        real(8), intent(in) :: alpha
        real(8), intent(inout) :: a(N), b(N), c(N)
        real(8), intent(out) :: phi(N)

        real(8) :: X_interface
        real(8), dimension(2) :: Dif = [1.0/3.0d0, 2.0/3.0d0]
        real(8), dimension(2) :: Sigma_a = [0.01, 0.02]
        real(8), dimension(2) :: S0 = [1.0, 1.0]
        real(8), dimension(2) :: dx
        real(8), dimension(2) :: L
        real(8), dimension(N) :: d, dxVis, Sigma_aVis
        integer :: ii, IDX

        X_interface = ((real(N_interface)/real(N)))
        dx = [X_interface/(N_interface-1), (1-X_interface)/(N-N_interface-1)]
        L = [sqrt(Dif(1)/Sigma_a(1)), sqrt(Dif(2)/Sigma_a(2))]

        do ii = 1, N
            if (ii <= N_interface) then
                IDX=1
            else
                IDX=2
            end if
            if (ii == N) then
                a(ii) = -2*Dif(IDX)/dx(IDX)**2
                b(ii) = Sigma_a(IDX) + 2.0*Dif(IDX)/dx(IDX)**2 + (1/dx(IDX))*(1-alpha)/(1+alpha)
                c(ii) = 0.0
                d(ii) = S0(IDX)
            elseif (ii == 1) then
                a(ii) = 0.0
                b(ii) = Sigma_a(IDX) + 2*Dif(IDX)/dx(IDX)**2 + (1/dx(IDX))*(1-alpha)/(1+alpha)
                c(ii) = -2*Dif(IDX)/dx(IDX)**2
                d(ii) = S0(IDX)
            else
                a(ii) = -Dif(IDX)/dx(IDX)**2
                b(ii) = Sigma_a(IDX) + 2.0*Dif(IDX)/dx(IDX)**2
                c(ii) = -Dif(IDX)/dx(IDX)**2
                d(ii) = S0(IDX)
            end if
            dxVis(ii) = dx(IDX)*ii
            Sigma_aVis(ii) = Sigma_a(IDX) 
        end do

        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

        ! Write all columns at once after the loop
        open(unit=991, file="flux.dat", status="replace", action="write")


        do ii = 1, N
            write(991,'(F10.5,2(1X,E15.8))') dxVis(ii), phi(ii)
        end do
        close(991)

    end subroutine build_matrix_A_heterogeneous_m1

        subroutine build_matrix_A_heterogeneous_m2(a, b, c, N, N_interface, x_interface, alpha, phi)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N, N_interface
        real(8), intent(in) :: alpha, x_interface
        real(8), intent(inout) :: a(N), b(N), c(N)
        real(8), intent(out) :: phi(N)

        real(8), dimension(N) :: Dif, Sigma_a, S0, L, d, phi_analytical, AA, x
        real(8), dimension(N-1) :: dx
        integer :: ii
        real(8) :: R_domain = 1

        Dif(1:N_interface-1) = 1.0/3.0d0 ;  Dif(N_interface:N) = 1.0/3.0d0
        Sigma_a(1:N_interface-1) = 0.01 ; Sigma_a(N_interface:N) = 0.01
        S0(1:N_interface-1) = 1.0 ; S0(N_interface:N) = 1.0

        ! First segment: 0 to x_interface
        dx(1:N_interface-1) = x_interface / real(N_interface-1, 8)

        ! Second segment: x_interface to 1.0
        dx(N_interface:N-1) = (1.0d0 - x_interface) / real(N - N_interface, 8)

        L = sqrt(Dif(1:N)/Sigma_a(1:N))

        x(1) = 0.0d0
        do ii = 2, N
            x(ii) = x(ii-1) + dx(ii-1)
        end do

        do ii = 1, N
            if (ii == 1) then
                a(ii) = 0.0d0
                b(ii) = ((Dif(ii)+Dif(ii))/dx(ii)+(Dif(ii)+Dif(ii+1))/dx(ii))/(dx(ii)+dx(ii))+Sigma_a(ii) + 1/dx(ii)*(1-alpha)/(1+alpha)
                c(ii) = -(Dif(ii)+Dif(ii+1))/(dx(ii)*(dx(ii)))
                d(ii) = S0(ii)
            elseif (ii == N) then
                a(ii) = -(Dif(ii-1)+Dif(ii))/(dx(ii-1)*dx(ii-1))
                b(ii) = ((Dif(ii)+Dif(ii-1))/dx(ii-1)+(Dif(ii)+Dif(ii))/dx(ii-1))/(dx(ii-1)+dx(ii-1))+Sigma_a(ii) + 1/dx(ii-1)*(1-alpha)/(1+alpha)
                c(ii) = 0.0d0
                d(ii) = S0(ii)
            else
                a(ii) = -(Dif(ii-1)+Dif(ii))/(dx(ii-1)*(dx(ii)+dx(ii-1)))
                b(ii) = ((Dif(ii)+Dif(ii-1))/dx(ii-1)+(Dif(ii)+Dif(ii+1))/dx(ii))/(dx(ii)+dx(ii-1))+Sigma_a(ii) 
                c(ii) = -(Dif(ii)+Dif(ii+1))/(dx(ii)*(dx(ii)+dx(ii-1)))
                d(ii) = S0(ii)
        end if

        end do

        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

        ! Write all columns at once after the loop
        open(unit=991, file="flux.dat", status="replace", action="write")
        do ii = 1, N
            AA(ii) = - ( ( S0(ii) / (2.0d0 * Sigma_a(ii)) ) * ( ( L(ii) / (2.0d0 * Dif(ii)) ) * sinh(R_domain / L(ii)) + 1.0d0 + cosh(R_domain / L(ii)) ) ) &
          / ( cosh(R_domain / L(ii)) + ( ( L(ii) / (4.0d0 * Dif(ii)) ) + ( Dif(ii) / L(ii) ) ) * sinh(R_domain / L(ii)) )
            !phi_analytical(ii) = AA(ii) * cosh( x(ii) / L(ii) ) + ( L(ii) / (2.0d0 * Dif(ii)) ) * ( AA(ii) + S0(ii) / Sigma_a(ii) ) * sinh( x(ii) / L(ii) ) + S0(ii) / Sigma_a(ii)
            write(991,'(F10.5,2(1X,E15.8))') x(ii), phi(ii), phi_analytical(ii)
        end do
        close(991)

    end subroutine build_matrix_A_heterogeneous_m2

end module m_Vacuum_BC_solver
