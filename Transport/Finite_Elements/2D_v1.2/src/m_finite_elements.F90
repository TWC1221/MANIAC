module m_finite_elements
    use m_constants
    use m_quadrature
    use m_types
    use m_constants, only: dp_EPSILON, check_nan_scalar, check_nan_array, check_nan_matrix
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

        if (allocated(FE%knots)) deallocate(FE%knots)
        allocate(FE%knots(2 * (FE%order + 1)))
        FE%knots = [ ( -1.0_dp, i=1, FE%order+1 ), ( 1.0_dp, i=1, FE%order+1 ) ]

        if (allocated(FE%basis_at_quad))        deallocate(FE%basis_at_quad)
        if (allocated(FE%dbasis_dxi))           deallocate(FE%dbasis_dxi)
        if (allocated(FE%dbasis_deta))          deallocate(FE%dbasis_deta)
        if (allocated(FE%basis_at_bound_quad))  deallocate(FE%basis_at_bound_quad)
        if (allocated(FE%dbasis_at_bound_quad)) deallocate(FE%dbasis_at_bound_quad)
        if (allocated(FE%face_node_map))        deallocate(FE%face_node_map)

        allocate(FE%basis_at_quad(Quad%n_points, FE%n_basis))
        allocate(FE%dbasis_dxi(Quad%n_points, FE%n_basis))
        allocate(FE%dbasis_deta(Quad%n_points, FE%n_basis))
        
        allocate(FE%basis_at_bound_quad(QuadBound%n_points, FE%order + 1))
        allocate(FE%dbasis_at_bound_quad(QuadBound%n_points, FE%order + 1))
        
        allocate(FE%face_node_map(FE%n_nodes_per_face, 4))

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
                    
                    FE%basis_at_quad(q, n) = FE_basis_1D(FE%knots, FE%order, i, Quad%xi(q)) * FE_basis_1D(FE%knots, FE%order, j, Quad%eta(q))

                    FE%dbasis_dxi(q, n) = FE_basis_derivative_1D(FE%knots, FE%order, i, Quad%xi(q)) * &
                                             FE_basis_1D(FE%knots, FE%order, j, Quad%eta(q))
                    
                    FE%dbasis_deta(q, n) = FE_basis_1D(FE%knots, FE%order, i, Quad%xi(q)) * &
                                              FE_basis_derivative_1D(FE%knots, FE%order, j, Quad%eta(q))
                    n = n + 1
                end do
            end do
        end do

        do q = 1, QuadBound%n_points
            do i = 1, FE%order + 1
                FE%basis_at_bound_quad(q, i) = FE_basis_1D(FE%knots, FE%order, i, QuadBound%xi(q))
                FE%dbasis_at_bound_quad(q, i) = FE_basis_derivative_1D(FE%knots, FE%order, i, QuadBound%xi(q))
            end do
        end do
        
    end subroutine InitialiseFiniteElements

    subroutine GetMapping(FE, q, elem_coords, dN_dx, dN_dy, detJ, elem_weights, basis_rational)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: q            
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: detJ
        real(dp), intent(in),  optional :: elem_weights(:)
        real(dp), intent(out), optional :: basis_rational(:)

        real(dp) :: J(2,2), invJ(2,2), dNdXiEta(2, FE%n_basis), N(FE%n_basis)
        real(dp) :: w_sum, dw_dxi, dw_deta

        N = FE%basis_at_quad(q, :)

        if (present(elem_weights)) then
            w_sum   = dot_product(N, elem_weights)
            dw_dxi  = dot_product(FE%dbasis_dxi(q, :), elem_weights)
            dw_deta = dot_product(FE%dbasis_deta(q, :), elem_weights)

            call check_nan_scalar(w_sum, "w_sum", "GetMapping, QP "//int_to_str(q))
            call check_nan_scalar(dw_dxi, "dw_dxi", "GetMapping, QP "//int_to_str(q))
            call check_nan_scalar(dw_deta, "dw_deta", "GetMapping, QP "//int_to_str(q))
            ! Quotient Rule for rational basis derivatives:
            dNdXiEta(1, :) = (FE%dbasis_dxi(q, :) * elem_weights * w_sum - (N * elem_weights) * dw_dxi) / (w_sum**2)
            dNdXiEta(2, :) = (FE%dbasis_deta(q, :) * elem_weights * w_sum - (N * elem_weights) * dw_deta) / (w_sum**2)
            
            if (present(basis_rational)) basis_rational = (N * elem_weights) / w_sum
        else
            dNdXiEta(1, :) = FE%dbasis_dxi(q, :)
            dNdXiEta(2, :) = FE%dbasis_deta(q, :)
            call check_nan_matrix(dNdXiEta, "dNdXiEta (non-rational)", "GetMapping, QP "//int_to_str(q))
            if (present(basis_rational)) basis_rational = N
        end if
        call check_nan_matrix(dNdXiEta, "dNdXiEta", "GetMapping, QP "//int_to_str(q))
        call check_nan_array(basis_rational, "basis_rational", "GetMapping, QP "//int_to_str(q))

        J = matmul(dNdXiEta, elem_coords)

        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)

        if (abs(detJ) < 1.0d-15) then
            dN_dx = 0.0_dp; dN_dy = 0.0_dp; detJ = 0.0_dp
            return
        end if

        invJ(1,1) =  J(2,2) / detJ
        invJ(1,2) = -J(1,2) / detJ
        invJ(2,1) = -J(2,1) / detJ
        invJ(2,2) =  J(1,1) / detJ

        dN_dx = invJ(1,1) * dNdXiEta(1, :) + invJ(2,1) * dNdXiEta(2, :)
        dN_dy = invJ(1,2) * dNdXiEta(1, :) + invJ(2,2) * dNdXiEta(2, :)

    end subroutine GetMapping

    subroutine GetArbitraryBasis(FE, xi, eta, N, elem_weights)
        type(t_finite), intent(in) :: FE
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: N(:)
        real(dp), intent(in), optional :: elem_weights(:)
        real(dp) :: w_sum
        integer :: i, j, nn
        nn = 1
        do j = 1, FE%order + 1
            do i = 1, FE%order + 1
                N(nn) = FE_basis_1D(FE%knots, FE%order, i, xi) * FE_basis_1D(FE%knots, FE%order, j, eta)
                nn = nn + 1
            end do
        end do

        if (present(elem_weights)) then
            N = N * elem_weights
            w_sum = sum(N)
            if (abs(w_sum) > dp_EPSILON) N = N / w_sum
        end if
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

    recursive function FE_basis_1D(knots, p, i, Xi) result(SF)
        real(dp), intent(in)    :: knots(:)
        integer, intent(in)     :: p
        integer, intent(in)     :: i
        real(dp), intent(in)    :: Xi
        real(dp)                :: SF

        SF = bspline_recurrence(knots, i, p, Xi)
    end function FE_basis_1D

    recursive function FE_basis_derivative_1D(knots, p, i, Xi) result(dSF)
        real(dp), intent(in)    :: knots(:)
        integer, intent(in)     :: p
        integer, intent(in)     :: i
        real(dp), intent(in)    :: Xi
        real(dp)                :: dSF
        real(dp)                :: denom1, denom2

        dSF = 0.0_dp
        denom1 = knots(i+p) - knots(i)
        denom2 = knots(i+p+1) - knots(i+1)

        if (denom1 > dp_EPSILON) dSF = dSF + p / denom1 * bspline_recurrence(knots, i, p-1, Xi)
        if (denom2 > dp_EPSILON) dSF = dSF - p / denom2 * bspline_recurrence(knots, i+1, p-1, Xi)
    end function FE_basis_derivative_1D

    recursive function bspline_recurrence(knots, i, p, xi) result(Ni)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: i, p
        real(dp), intent(in) :: xi
        real(dp) :: Ni, d1, d2

        if (p == 0) then
            call check_nan_scalar(knots(i), "knots(i)", "bspline_recurrence")
            call check_nan_scalar(knots(i+1), "knots(i+1)", "bspline_recurrence")
            if (xi >= knots(i) .and. xi < knots(i+1)) then
                Ni = 1.0_dp
            else if (xi == knots(i+1) .and. i+1 == size(knots)) then
                Ni = 1.0_dp
            else
                Ni = 0.0_dp
            end if
        else
            Ni = 0.0_dp
            d1 = knots(i+p) - knots(i)
            d2 = knots(i+p+1) - knots(i+1)
            call check_nan_scalar(d1, "d1", "bspline_recurrence")
            call check_nan_scalar(d2, "d2", "bspline_recurrence")
            if (d1 > dp_EPSILON) Ni = Ni + (xi - knots(i)) / d1 * bspline_recurrence(knots, i, p-1, xi)
            if (d2 > dp_EPSILON) Ni = Ni + (knots(i+p+1) - xi) / d2 * bspline_recurrence(knots, i+1, p-1, xi)
            call check_nan_scalar(Ni, "Ni (recursive)", "bspline_recurrence")
        end if
    end function bspline_recurrence

    function outer_product(a, b) result(mat)
        real(dp), intent(in)  :: a(:), b(:)
        real(dp)              :: mat(size(a,1), size(b,1))
        mat = spread(a, dim=2, ncopies=size(b,1)) * spread(b, dim=1, ncopies=size(a,1))
    end function outer_product
    
end module m_finite_elements