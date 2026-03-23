module m_finite_elements
    use m_constants
    use m_quadrature
    use m_types
    implicit none
    public
        
    contains

    subroutine InitialiseFiniteElements(FE, Quad, QuadBound, write_nodal_map)
        type(t_finite), intent(inout) :: FE
        type(t_quadrature), intent(in) :: Quad, QuadBound
        integer :: q, i, j, k, n, idx, Poly
        logical :: write_nodal_map

        Poly = FE%order
        FE%n_basis = (Poly + 1)**3
        
        call finite_1D_positions(FE, Poly)

        if (allocated(FE%N))       deallocate(FE%N)
        if (allocated(FE%dN_dxi))  deallocate(FE%dN_dxi)
        if (allocated(FE%dN_deta)) deallocate(FE%dN_deta)
        if (allocated(FE%dN_dzeta)) deallocate(FE%dN_dzeta)
        if (allocated(FE%N_mat))   deallocate(FE%N_mat)
        if (allocated(FE%N_B))     deallocate(FE%N_B)
        if (allocated(FE%dN_B_xi)) deallocate(FE%dN_B_xi)
        if (allocated(FE%dN_B_eta)) deallocate(FE%dN_B_eta)
        if (allocated(FE%p))       deallocate(FE%p)

        allocate(FE%N(Quad%NoPoints, FE%n_basis))
        allocate(FE%dN_dxi(Quad%NoPoints, FE%n_basis))
        allocate(FE%dN_deta(Quad%NoPoints, FE%n_basis))
        allocate(FE%dN_dzeta(Quad%NoPoints, FE%n_basis))
        allocate(FE%N_mat(Quad%NoPoints, FE%n_basis, FE%n_basis))
        
        allocate(FE%N_B(QuadBound%NoPoints, (Poly + 1)**2))
        allocate(FE%dN_B_xi(QuadBound%NoPoints, (Poly + 1)**2))
        allocate(FE%dN_B_eta(QuadBound%NoPoints, (Poly + 1)**2))
        allocate(FE%p(FE%n_basis))


        call GetVTKOrderVector(FE)
        if (write_nodal_map) then
        print*, "ELement Nodal INDEX Map:"
        print*, "---------------------------------"
        do j = Poly + 1, 1, -1
            write(*, '(A2, I2, A3)', advance='no') "j=", j, " |"
            do i = 1, Poly + 1
                idx = i + (j-1)*(Poly + 1)
                write(*, '(I4)', advance='no') FE%p(idx)
            end do
            write(*,*)
        end do
        print*, "---------------------------------"
        write(*, '(A, I2)') " Element Order: ", Poly
        end if

        do q = 1, Quad%NoPoints
            n = 1
            do k = 1, Poly + 1         ! z-direction
                do j = 1, Poly + 1     ! y-direction
                    do i = 1, Poly + 1 ! x-direction
                        
                        FE%N(q, FE%p(n)) = FE_basis_1D(FE, Poly, i, Quad%Xi(q)) * &
                                           FE_basis_1D(FE, Poly, j, Quad%Eta(q)) * &
                                           FE_basis_1D(FE, Poly, k, Quad%Zeta(q))

                        FE%dN_dxi(q, FE%p(n)) = FE_basis_derivative_1D(FE, Poly, i, Quad%Xi(q)) * &
                                                FE_basis_1D(FE, Poly, j, Quad%Eta(q)) * &
                                                FE_basis_1D(FE, Poly, k, Quad%Zeta(q))
                        
                        FE%dN_deta(q, FE%p(n)) = FE_basis_1D(FE, Poly, i, Quad%Xi(q)) * &
                                                FE_basis_derivative_1D(FE, Poly, j, Quad%Eta(q)) * &
                                                FE_basis_1D(FE, Poly, k, Quad%Zeta(q))

                        FE%dN_dzeta(q, FE%p(n)) = FE_basis_1D(FE, Poly, i, Quad%Xi(q)) * &
                                                  FE_basis_1D(FE, Poly, j, Quad%Eta(q)) * &
                                                  FE_basis_derivative_1D(FE, Poly, k, Quad%Zeta(q))
                        n = n + 1
                    end do
                end do
            end do
            FE%N_mat(q,:,:) = outer_product(FE%N(q,:), FE%N(q,:))
        end do

        do q = 1, QuadBound%NoPoints
            n = 1
            do j = 1, Poly + 1
                do i = 1, Poly + 1
                    FE%N_B(q, n) = FE_basis_1D(FE, Poly, i, QuadBound%Xi(q)) * &
                                   FE_basis_1D(FE, Poly, j, QuadBound%Eta(q))
                    FE%dN_B_xi(q, n) = FE_basis_derivative_1D(FE, Poly, i, QuadBound%Xi(q)) * &
                                       FE_basis_1D(FE, Poly, j, QuadBound%Eta(q))
                    FE%dN_B_eta(q, n) = FE_basis_1D(FE, Poly, i, QuadBound%Xi(q)) * &
                                        FE_basis_derivative_1D(FE, Poly, j, QuadBound%Eta(q))
                    n = n + 1
                end do
            end do
        end do
    end subroutine InitialiseFiniteElements

    subroutine GetMapping(FE, q, elem_coords, dN_dx, dN_dy, dN_dz, detJ)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: q            
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:), dN_dz(:)
        real(dp), intent(out)        :: detJ

        real(dp) :: J(3,3), invJ(3,3)
        integer :: i

        J = 0.0_dp
        do i = 1, FE%n_basis
            J(1,1) = J(1,1) + FE%dN_dxi(q, i)  * elem_coords(i, 1)
            J(1,2) = J(1,2) + FE%dN_dxi(q, i)  * elem_coords(i, 2)
            J(1,3) = J(1,3) + FE%dN_dxi(q, i)  * elem_coords(i, 3)
            J(2,1) = J(2,1) + FE%dN_deta(q, i) * elem_coords(i, 1)
            J(2,2) = J(2,2) + FE%dN_deta(q, i) * elem_coords(i, 2)
            J(2,3) = J(2,3) + FE%dN_deta(q, i) * elem_coords(i, 3)
            J(3,1) = J(3,1) + FE%dN_dzeta(q, i) * elem_coords(i, 1)
            J(3,2) = J(3,2) + FE%dN_dzeta(q, i) * elem_coords(i, 2)
            J(3,3) = J(3,3) + FE%dN_dzeta(q, i) * elem_coords(i, 3)
        end do

        detJ = J(1,1)*(J(2,2)*J(3,3)-J(3,2)*J(2,3)) - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))

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

        invJ(1,1) = (J(2,2)*J(3,3) - J(2,3)*J(3,2)) / detJ
        invJ(1,2) = (J(1,3)*J(3,2) - J(1,2)*J(3,3)) / detJ
        invJ(1,3) = (J(1,2)*J(2,3) - J(1,3)*J(2,2)) / detJ

        invJ(2,1) = (J(2,3)*J(3,1) - J(2,1)*J(3,3)) / detJ
        invJ(2,2) = (J(1,1)*J(3,3) - J(1,3)*J(3,1)) / detJ
        invJ(2,3) = (J(1,3)*J(2,1) - J(1,1)*J(2,3)) / detJ

        invJ(3,1) = (J(2,1)*J(3,2) - J(2,2)*J(3,1)) / detJ
        invJ(3,2) = (J(1,2)*J(3,1) - J(1,1)*J(3,2)) / detJ
        invJ(3,3) = (J(1,1)*J(2,2) - J(1,2)*J(2,1)) / detJ

        do i = 1, FE%n_basis
            dN_dx(i) = invJ(1,1) * FE%dN_dxi(q, i) + invJ(1,2) * FE%dN_deta(q, i) + invJ(1,3) * FE%dN_dzeta(q, i)
            dN_dy(i) = invJ(2,1) * FE%dN_dxi(q, i) + invJ(2,2) * FE%dN_deta(q, i) + invJ(2,3) * FE%dN_dzeta(q, i)
            dN_dz(i) = invJ(3,1) * FE%dN_dxi(q, i) + invJ(3,2) * FE%dN_deta(q, i) + invJ(3,3) * FE%dN_dzeta(q, i)
        end do

    end subroutine GetMapping

subroutine GetVTKOrderVector(FE)
        type(t_finite), intent(inout) :: FE
        integer :: P, n, i, j, idx, L, curP
        integer, allocatable :: temp_map(:,:)
        
        P = FE%order
        if (allocated(FE%p)) deallocate(FE%p)
        allocate(FE%p(FE%n_basis))
        allocate(temp_map(P+1, P+1))

        n = 1
        L = 1
        curP = P
        
        if (P == 1) then
            FE%p = [1, 2, 4, 3, 5, 6, 8, 7]
            deallocate(temp_map)
            return
        end if

        do i = 1, FE%n_basis; FE%p(i) = i; end do

        do while (curP >= 1)
            temp_map(L, L)             = n; n = n + 1 ! BL
            temp_map(L + curP, L)      = n; n = n + 1 ! BR
            temp_map(L + curP, L + curP)= n; n = n + 1 ! TR
            temp_map(L, L + curP)      = n; n = n + 1 ! TL

            if (curP > 1) then
                do i = L + 1, L + curP - 1
                    temp_map(i, L) = n; n = n + 1 ! Bottom
                end do
                do j = L + 1, L + curP - 1
                    temp_map(L + curP, j) = n; n = n + 1 ! Right
                end do
                do i = L + curP - 1, L + 1, -1
                    temp_map(i, L + curP) = n; n = n + 1 ! Top
                end do
                do j = L + curP - 1, L + 1, -1
                    temp_map(L, j) = n; n = n + 1 ! Left
                end do
            end if

            L = L + 1
            curP = curP - 2
        end do

        if (curP == 0) temp_map(L, L) = n

        idx = 1
        do j = 1, P+1
            do i = 1, P+1
                FE%p(idx) = temp_map(i, j)
                idx = idx + 1
            end do
        end do
        deallocate(temp_map)
    end subroutine GetVTKOrderVector

    subroutine finite_1D_positions(FE, PolyOrder)
        type(t_finite), intent(inout) :: FE
        integer, intent(in)         :: PolyOrder
        integer :: i
        real(dp) :: dxi

        if (allocated(FE%Xi)) deallocate(FE%Xi)
        allocate(FE%Xi(PolyOrder + 1))

        if (PolyOrder == 1) then
            FE%Xi(1) = -1.0_dp
            FE%Xi(2) =  1.0_dp
            return
        end if

        dxi = 2.0_dp / real(PolyOrder, dp)
        do i = 1, PolyOrder + 1
            FE%Xi(i) = -1.0_dp + (i-1)*dxi
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