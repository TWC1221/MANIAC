module m_finite_elements
    use m_constants
    use m_quadrature
    use m_types
    implicit none
    public
        
    contains

    subroutine InitFE(FE, Quad, QuadBound)
        type(t_finite), intent(inout) :: FE
        type(t_quadrature), intent(in) :: Quad, QuadBound
        integer :: gp, i, j, n, Poly
        
        Poly = FE%order
        FE%n_basis = (Poly + 1)**2
        
        call finite_1D_positions(FE, Poly)

        if (allocated(FE%N))      deallocate(FE%N)
        if (allocated(FE%dN_dxi))  deallocate(FE%dN_dxi)
        if (allocated(FE%dN_deta)) deallocate(FE%dN_deta)
        if (allocated(FE%N_mat))   deallocate(FE%N_mat)
        if (allocated(FE%N_B))     deallocate(FE%N_B)
        if (allocated(FE%dN_B))    deallocate(FE%dN_B)
        if (allocated(FE%p)) deallocate(FE%p)

        allocate(FE%N(Quad%NoPoints, FE%n_basis))
        allocate(FE%dN_dxi(Quad%NoPoints, FE%n_basis))
        allocate(FE%dN_deta(Quad%NoPoints, FE%n_basis))
        allocate(FE%N_mat(Quad%NoPoints, FE%n_basis, FE%n_basis))
        
        allocate(FE%N_B(QuadBound%NoPoints, Poly + 1))
        allocate(FE%dN_B(QuadBound%NoPoints, Poly + 1))
        allocate(FE%p(FE%n_basis))

        select case (FE%n_basis)
        case (4)  ; FE%p = [1, 2, 4, 3]                  ! Linear
        case (9)  ; FE%p = [1, 5, 2, 8, 9, 6, 4, 7, 3]   ! Quadratic
        case (16) ; FE%p = [1, 5, 6, 2, 12, 13, 14, 7, 11, 16, 15, 8, 4, 10, 9, 3] ! Cubic
        end select
            
        do gp = 1, Quad%NoPoints
            n = 1
            do j = 1, Poly + 1     ! y-direction
                do i = 1, Poly + 1 ! x-direction
                    
                    FE%N(gp, FE%p(n)) = FE_basis_1D(FE, Poly, i, Quad%Xi(gp)) * &
                                        FE_basis_1D(FE, Poly, j, Quad%Eta(gp))
                    
                    FE%dN_dxi(gp, FE%p(n)) = FE_basis_derivative_1D(FE, Poly, i, Quad%Xi(gp)) * &
                                             FE_basis_1D(FE, Poly, j, Quad%Eta(gp))
                    
                    FE%dN_deta(gp, FE%p(n)) = FE_basis_1D(FE, Poly, i, Quad%Xi(gp)) * &
                                              FE_basis_derivative_1D(FE, Poly, j, Quad%Eta(gp))
                    n = n + 1
                end do
            end do
            FE%N_mat(gp,:,:) = outer_product(FE%N(gp,:), FE%N(gp,:))
        end do

        do gp = 1, QuadBound%NoPoints
            do i = 1, Poly + 1
                FE%N_B(gp, i) = FE_basis_1D(FE, Poly, i, QuadBound%Xi(gp))
                FE%dN_B(gp, i) = FE_basis_derivative_1D(FE, Poly, i, QuadBound%Xi(gp))
            end do
        end do
    end subroutine InitFE

    subroutine GetMapping(FE, gp, elem_coords, dN_dx, dN_dy, detJ)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: gp            
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: detJ

        real(dp) :: J(2,2), invJ(2,2)
        integer :: i

        J = 0.0_dp
        do i = 1, FE%n_basis
            J(1,1) = J(1,1) + FE%dN_dxi(gp, i)  * elem_coords(i, 1)
            J(1,2) = J(1,2) + FE%dN_dxi(gp, i)  * elem_coords(i, 2)
            J(2,1) = J(2,1) + FE%dN_deta(gp, i) * elem_coords(i, 1)
            J(2,2) = J(2,2) + FE%dN_deta(gp, i) * elem_coords(i, 2)
        end do

        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)

        if (detJ <= 1.0d-10) then
            print *, " "
            print *, ">>> GEOMETRY ERROR: detJ =", detJ
            print *, ">>> Check node ordering for element coordinates below:"
            do i = 1, FE%n_basis
                write(*, '(A, I2, A, F10.5, A, F10.5, A)') "Node ", i, ": (", &
                    elem_coords(i,1), ", ", elem_coords(i,2), ")"
            end do
            print *, ">>> Check if Nodes 2 and 4 are swapped in your map!"
            print *, " "
            stop "FATAL: Mapping failed due to non-positive Jacobian."
        end if

        invJ(1,1) =  J(2,2) / detJ
        invJ(1,2) = -J(1,2) / detJ
        invJ(2,1) = -J(2,1) / detJ
        invJ(2,2) =  J(1,1) / detJ

        do i = 1, FE%n_basis
            dN_dx(i) = invJ(1,1) * FE%dN_dxi(gp, i) + invJ(1,2) * FE%dN_deta(gp, i)
            dN_dy(i) = invJ(2,1) * FE%dN_dxi(gp, i) + invJ(2,2) * FE%dN_deta(gp, i)
        end do

    end subroutine GetMapping

    subroutine finite_1D_positions(Finite, PolyOrder)
        type(t_finite), intent(inout) :: Finite
        integer, intent(in)         :: PolyOrder
        integer :: i
        real(dp) :: dxi

        if (allocated(Finite%Xi)) deallocate(Finite%Xi)
        allocate(Finite%Xi(PolyOrder + 1))

        if (PolyOrder == 1) then
            Finite%Xi(1) = -1.0_dp
            Finite%Xi(2) =  1.0_dp
            return
        end if

        dxi = 2.0_dp / real(PolyOrder, dp)
        do i = 1, PolyOrder + 1
            Finite%Xi(i) = -1.0_dp + (i-1)*dxi
        end do
    end subroutine finite_1D_positions

    real(dp) function FE_basis_1D(FE, Poly, i, Xi) result(SF)
        type(t_finite)          :: FE
        integer, intent(in)     :: i, Poly
        real(dp), intent(in)    :: Xi
        integer                 :: p
        SF = 1.0_dp
        do p = 1, 1 + Poly
            if (p /= i) SF = SF * (Xi - FE%Xi(p)) / (FE%Xi(i) - FE%Xi(p))
        end do
    end function FE_basis_1D

    real(dp) function FE_basis_derivative_1D(FE, Poly, i, Xi) result(dSF)
        type(t_finite)          :: FE
        integer, intent(in)     :: i, Poly
        real(dp), intent(in)    :: Xi
        real(dp)                :: TempValue
        integer                 :: p, q
        dSF = 0.0_dp
        do p = 1, 1 + Poly
            if (p /= i) then
                TempValue = 1.0_dp / (FE%Xi(i) - FE%Xi(p))
                do q = 1, 1 + Poly
                    if (q /= i .and. q /= p) then
                        TempValue = TempValue * (Xi - FE%Xi(q)) / (FE%Xi(i) - FE%Xi(q))
                    end if
                end do
                dSF = dSF + TempValue
            end if
        end do
    end function FE_basis_derivative_1D

    function outer_product(a, b) result(mat)
        real(8), intent(in)  :: a(:), b(:)
        real(8)              :: mat(size(a,1), size(b,1))
        mat = spread(a, dim=2, ncopies=size(b,1)) * spread(b, dim=1, ncopies=size(a,1))
    end function outer_product
end module m_finite_elements