module m_finite_elements
    use m_constants
    use m_quadrature
    use m_types
    implicit none
    public
        
    contains
    
    subroutine InitialiseFiniteElements(FE, QuadBound, Quad)
        type(t_finite), intent(inout) :: FE
        type(t_quadrature), intent(in) :: Quad, QuadBound
        integer :: q, i, j, n
        
        FE%n_nodes_per_face = FE%order + 1
        FE%n_basis = (FE%order + 1)**2

        call finite_1D_positions(FE)

        if (allocated(FE%basis_at_quad))           deallocate(FE%basis_at_quad)
        if (allocated(FE%dbasis_dxi))              deallocate(FE%dbasis_dxi)
        if (allocated(FE%dbasis_deta))             deallocate(FE%dbasis_deta)
        if (allocated(FE%basis_at_bound_quad))     deallocate(FE%basis_at_bound_quad)
        if (allocated(FE%dbasis_at_bound_quad))    deallocate(FE%dbasis_at_bound_quad)
        if (allocated(FE%face_node_map))           deallocate(FE%face_node_map)

        allocate(FE%basis_at_quad(Quad%n_points, FE%n_basis))
        allocate(FE%dbasis_dxi(Quad%n_points, FE%n_basis))
        allocate(FE%dbasis_deta(Quad%n_points, FE%n_basis))
        allocate(FE%basis_at_bound_quad(QuadBound%n_points, FE%order + 1))
        allocate(FE%dbasis_at_bound_quad(QuadBound%n_points, FE%order + 1))
        allocate(FE%face_node_map(FE%n_nodes_per_face, 4))

        ! Local indices for a tensor-product element with (p+1)x(p+1) nodes
        do i = 1, FE%n_nodes_per_face
            FE%face_node_map(i, 1) = i
            FE%face_node_map(i, 2) = (i-1)*(FE%order+1) + (FE%order+1)
            FE%face_node_map(i, 3) = FE%n_basis - (i-1)
            FE%face_node_map(i, 4) = (FE%n_nodes_per_face - i)*(FE%order+1) + 1
        end do

        do q = 1, Quad%n_points
            n = 1
            do j = 1, FE%order + 1     
                do i = 1, FE%order + 1 
                    
                    FE%basis_at_quad(q, n) = FE_basis_1D(FE, i, Quad%xi(q)) * &
                                        FE_basis_1D(FE, j, Quad%eta(q))

                    FE%dbasis_dxi(q, n) = FE_basis_derivative_1D(FE, i, Quad%xi(q)) * &
                                             FE_basis_1D(FE, j, Quad%eta(q))
                    
                    FE%dbasis_deta(q, n) = FE_basis_1D(FE, i, Quad%xi(q)) * &
                                              FE_basis_derivative_1D(FE, j, Quad%eta(q))
                    n = n + 1
                end do
            end do
        end do

        do q = 1, QuadBound%n_points
            do i = 1, FE%order + 1
                FE%basis_at_bound_quad(q, i) = FE_basis_1D(FE, i, QuadBound%xi(q))
                FE%dbasis_at_bound_quad(q, i) = FE_basis_derivative_1D(FE, i, QuadBound%xi(q))
            end do
        end do
    end subroutine InitialiseFiniteElements

    subroutine GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: q            
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: detJ

        real(dp) :: J(2,2), invJ(2,2), dNdXiEta(2, FE%n_basis)
        integer  :: i ! For error message loop

        dNdXiEta(1, :) = FE%dbasis_dxi(q, :)
        dNdXiEta(2, :) = FE%dbasis_deta(q, :)
        J = matmul(dNdXiEta, elem_coords)

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

        dN_dx = invJ(1,1) * FE%dbasis_dxi(q, :) + invJ(1,2) * FE%dbasis_deta(q, :)
        dN_dy = invJ(2,1) * FE%dbasis_dxi(q, :) + invJ(2,2) * FE%dbasis_deta(q, :)

    end subroutine GetMapping

    subroutine GetArbitraryBasis(FE, xi, eta, N)
        type(t_finite), intent(in) :: FE
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: N(:)
        integer :: i, j, nn
        nn = 1
        do j = 1, FE%order + 1
            do i = 1, FE%order + 1
                N(nn) = FE_basis_1D(FE, i, xi) * FE_basis_1D(FE, j, eta)
                nn = nn + 1
            end do
        end do
    end subroutine GetArbitraryBasis

    subroutine finite_1D_positions(FE)
        type(t_finite), intent(inout) :: FE
        integer :: i
        real(dp) :: dxi

        if (allocated(FE%node_roots)) deallocate(FE%node_roots)
        allocate(FE%node_roots(FE%order + 1))

        if (FE%order == 1) then
            FE%node_roots(1) = -1.0_dp
            FE%node_roots(2) =  1.0_dp
            return
        end if

        dxi = 2.0_dp / real(FE%order, dp)
        do i = 1, FE%order + 1
            FE%node_roots(i) = -1.0_dp + (i-1)*dxi
        end do
    end subroutine finite_1D_positions

    real(dp) function FE_basis_1D(FE, i, Xi) result(SF)
        type(t_finite)          :: FE
        integer, intent(in)     :: i
        real(dp), intent(in)    :: Xi
        integer                 :: p
        SF = 1.0_dp
        do p = 1, 1 + FE%order
            if (p /= i) SF = SF * (Xi - FE%node_roots(p)) / (FE%node_roots(i) - FE%node_roots(p))
        end do
    end function FE_basis_1D

    real(dp) function FE_basis_derivative_1D(FE, i, Xi) result(dSF)
        type(t_finite)          :: FE
        integer, intent(in)     :: i
        real(dp), intent(in)    :: Xi
        real(dp)                :: TempValue
        integer                 :: p, q
        dSF = 0.0_dp
        do p = 1, 1 + FE%order
            if (p /= i) then
                TempValue = 1.0_dp / (FE%node_roots(i) - FE%node_roots(p))
                do q = 1, 1 + FE%order
                    if (q /= i .and. q /= p) then
                        TempValue = TempValue * (Xi - FE%node_roots(q)) / (FE%node_roots(i) - FE%node_roots(q))
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