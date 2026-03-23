module m_utilib

!------------------------------------------------------------------------!
!! Date of creation:                                                    -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 28/08/2025    C. Jones         **************                         -!
!------------------------------------------------------------------------!
!
!! Module:
!   This module contains utility functions and subroutines for element
!   computations, including geometry, integration, and interpolation.
!
    use m_constants
    use m_quadrature
    implicit none
    public

contains


!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Matrix operations                                                    -!
!------------------------------------------------------------------------!  

subroutine inv_2x2_matrix(matrix, inverse, det)

    real(8), intent(in)     :: matrix(:,:)
    real(8), intent(out)    :: inverse(:,:)
    real(8)                 :: det
    
    ! Check if the determinant is nonzero (matrix is invertible)
    if (det == 0.0_dp) then
        write(*,*) 'Error: Matrix is singular. Cannot compute inverse.'
        stop
    end if

    ! Calculate inverse
    inverse(1,1) = matrix(2,2) / det
    inverse(1,2) = -matrix(1,2) / det
    inverse(2,1) = -matrix(2,1) / det
    inverse(2,2) = matrix(1,1) / det

end subroutine inv_2x2_matrix

subroutine inv_3x3_matrix(A, invA, det)

    real(dp), dimension(:,:)    :: A
    real(dp), dimension(:,:)    :: invA
    real(dp)                    :: det
    real(dp), dimension(3,3)    :: adj
    

    ! Calculate the adjugate matrix
    adj(1,1) =   (A(2,2) * A(3,3) - A(2,3) * A(3,2))
    adj(1,2) = - (A(1,2) * A(3,3) - A(1,3) * A(3,2))
    adj(1,3) =   (A(1,2) * A(2,3) - A(1,3) * A(2,2))
    
    adj(2,1) = - (A(2,1) * A(3,3) - A(2,3) * A(3,1))
    adj(2,2) =   (A(1,1) * A(3,3) - A(1,3) * A(3,1))
    adj(2,3) = - (A(1,1) * A(2,3) - A(1,3) * A(2,1))
    
    adj(3,1) =   (A(2,1) * A(3,2) - A(2,2) * A(3,1))
    adj(3,2) = - (A(1,1) * A(3,2) - A(1,2) * A(3,1))
    adj(3,3) =   (A(1,1) * A(2,2) - A(1,2) * A(2,1))

    ! Calculate the inverse of A
    invA = adj / det

end subroutine inv_3x3_matrix

real(dp) function calc_3x3_det(A) result(det)
    real(dp), dimension(:,:)    :: A  ! Input matrix

    ! Calculate the determinant using the rule of Sarrus for 3x3 matrices
    det = A(1,1) * (A(2,2) * A(3,3) - A(2,3) * A(3,2)) &
        - A(1,2) * (A(2,1) * A(3,3) - A(2,3) * A(3,1)) &
        + A(1,3) * (A(2,1) * A(3,2) - A(2,2) * A(3,1))

end function calc_3x3_det

function cross_product(a, b) result(c)

    real(kind=8), intent(in)    :: a(3), b(3)   ! Input vectors (3D)
    real(kind=8)                :: c(3)         ! Resulting cross product vector

    ! Compute the cross product
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

end function cross_product

function outer_product(a, b) result(mat)
    real(8), intent(in)  :: a(:), b(:)
    real(8)              :: mat(size(a,1), size(b,1))
    
    mat = spread(a, dim=2, ncopies=size(b,1)) * spread(b, dim=1, ncopies=size(a,1))
end function outer_product

function Invert_Matrix(A, N) result(AInv)

    integer, intent(in)                         :: N
    real(dp), dimension(N, N), intent(inout)    :: A
    real(dp), dimension(N, N)                   :: AInv

    integer                     :: INFO
    integer, dimension(N)       :: IPIV
    integer                     :: LWORK
    real(dp), dimension(N*N)    :: WORK
    real(dp), dimension(N,N)    :: TEMP

    LWORK = N * N
    TEMP = A

    ! Initialize INFO
    INFO = 0
    
    ! Call DGETRF to perform LU factorization
    call DGETRF(N, N, TEMP, N, IPIV, INFO)
    if (INFO /= 0) then
        print*, 'Error in LU factorization.'
        return
    endif

    ! Call DGETRI to compute inverse
    call DGETRI(N, TEMP, N, IPIV, WORK, LWORK, INFO)
    if (INFO /= 0) then
        print*, 'Error in computing inverse matrix.'
        return
    endif
    ! Copy the inverted matrix back to A
    AInv = TEMP
    
end function Invert_Matrix

subroutine print_matrix_real(mat)
    real(dp), dimension(:,:)    :: mat
    integer                     :: i, j

    do i = 1, size(mat,dim=1)
            write(*,'(A2)',advance="no") "[ "
            do j = 1,size(mat,dim=2)
                if (abs(mat(i,j)) < 1e-9) then
                    write(*,'(f10.0,A1)',advance="no") 0.0," "
                else
                    write(*,'(f10.5,A1)',advance="no") Mat(i,j)," "
                end if
            end do
            write(*,'(A1)') "]"
        end do
        write(*,*) ''

end subroutine print_matrix_real

subroutine print_matrix_int(mat)
    integer, dimension(:,:)     :: mat
    integer                     :: i, j

    do i = 1, size(mat,dim=1)
            write(*,'(A2)',advance="no") "[ "
            do j = 1,size(mat,dim=2)
                ! if (abs(mat(i,j)) == 0) then
                !     write(*,'(I2,A1)',advance="no") 0.0," "
                ! else
                    write(*,'(I2,A1)',advance="no") Mat(i,j)," "
                ! end if
            end do
            write(*,'(A1)') "]"
        end do
        write(*,*) ''
        
end subroutine print_matrix_int


!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Interpolation function                                               -!
!------------------------------------------------------------------------!  


REAL(kind =dp) function interpolate_func(x1,x2,x3,val1,val2)
    real(KIND = dp) :: x1,x2,x3,val1,val2
    interpolate_func = ((val2-val1)*(x3-x1))/(x2-x1)+val1
end function

REAL(kind =dp) function twoD_interpolate_func(x1,x2,x3,y1,y2,y3,val11,val12,val21,val22)
    real(KIND = dp) :: x1,x2,x3,y1,y2,y3,val11,val12,val21,val22,av1,av2,av3,av4
    av1 = interpolate_func(x1,x2,x3,val11,val12)
    av2 = interpolate_func(x1,x2,x3,val21,val22)
    av3 = interpolate_func(y1,y2,y3,val11,val21)
    av4 = interpolate_func(y1,y2,y3,val12,val22)
    twoD_interpolate_func = (av1+av2+av3+av4)/4
end function

!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Geometric functions                                                  -!
!------------------------------------------------------------------------! 

logical function is_anticlockwise(nodes)

    real(dp), intent(in)    :: nodes(2,2)
    real(dp)                :: cross

    ! Compute z-component of cross product between origin->node1 and origin->node2
    cross = nodes(1,1) * nodes(2,2) - nodes(2,1) * nodes(1,2)

    ! True if anticlockwise
    is_anticlockwise = (cross > 0.0_dp)

end function is_anticlockwise

end module m_utilib