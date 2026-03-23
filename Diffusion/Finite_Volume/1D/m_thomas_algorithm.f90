module m_thomas_algorithm
!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Contains the subroutines required to solve the tridiagonal matrix    -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 10/10/2025    T. Charlton      Original code                          -!
!------------------------------------------------------------------------! 

implicit none
contains
    !------------------------------------------------------------!
    ! Standard Thomas algorithm for tridiagonal systems         !
    !------------------------------------------------------------!

    subroutine thomas_algorithm(a, b, c, d, x, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: a(n-1), b(n), c(n-1), d(n)
        real(8), intent(out) :: x(n)
        real(8), allocatable :: c_star(:), d_star(:)
        integer :: ii

        allocate(c_star(n-1))
        allocate(d_star(n))

        ! Forward sweep
        c_star(1) = c(1)/b(1)
        d_star(1) = d(1)/b(1)

        do ii = 2, n-1
            c_star(ii) = c(ii)/(b(ii)-a(ii-1)*c_star(ii-1))
        end do

        do ii = 2, n
            d_star(ii) = (d(ii)-a(ii-1)*d_star(ii-1))/(b(ii)-a(ii-1)*c_star(ii-1))
        end do

        ! Back substitution
        x(n) = d_star(n)
        do ii = n-1, 1, -1
            x(ii) = d_star(ii) - c_star(ii)*x(ii+1)
        end do

        deallocate(c_star, d_star)

    end subroutine thomas_algorithm

end module