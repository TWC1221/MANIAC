module m_multigroup_1D_diffusion
implicit none
contains
    subroutine build_matrix_multigroup(N, alpha, G, phi_ptr, dx, Dif, A_scatter, S0, gg, Sigma_a)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N, G, gg
        real(8), intent(in) :: alpha
        real(8), intent(inout) :: phi_ptr(1:N), dx(N-1), Dif, A_scatter(1:G,1:G), S0(N), Sigma_a

        real(8), dimension(N) :: a, b, c, d, x
        integer :: ii

        x(1) = 0.0d0
        do ii = 2, N
            x(ii) =  sum(dx(1:ii-1))
        end do

        do ii = 1, N
            if (ii == 1) then
                a(ii) = 0.0d0
                b(ii) = (2*Dif)/(dx(ii)**2) + 1/dx(ii)*(1-alpha)/(1+alpha) + sum(A_scatter(gg,1:G)) - A_scatter(gg,gg) + Sigma_a !+ VLARGE_NUMBER
                c(ii) = -(2*Dif)/(dx(ii)**2)
                d(ii) = S0(ii)

            elseif (ii == N) then
                a(ii) = -(2*Dif)/(dx(ii-1)**2)
                b(ii) = (2*Dif)/(dx(ii-1)**2) + 1/dx(ii-1)*(1-alpha)/(1+alpha) + sum(A_scatter(gg,1:G)) - A_scatter(gg,gg) + Sigma_a !+ VLARGE_NUMBER
                c(ii) = 0.0d0
                d(ii) = S0(ii)
            else
                a(ii) = -(2*Dif)/(dx(ii-1)*(dx(ii)+dx(ii-1)))
                b(ii) = ((2*Dif)/dx(ii-1)+(2*Dif)/dx(ii))/(dx(ii)+dx(ii-1)) + sum(A_scatter(gg,1:G)) - A_scatter(gg,gg) + Sigma_a
                c(ii) = -(2*Dif)/(dx(ii)*(dx(ii)+dx(ii-1)))
                d(ii) = S0(ii)
            end if
        end do
        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi_ptr, N)

    end subroutine

end module m_multigroup_1D_diffusion