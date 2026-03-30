module m_finite_elements
    use m_constants
    use m_quadrature
    use m_types
    implicit none
    public
        
    contains
    
    subroutine InitialiseFiniteElements(FE, Quad2D, Quad3D)
        type(t_finite), intent(inout) :: FE
        type(t_quadrature), intent(in) :: Quad3D, Quad2D
        integer :: q, i, j, k, n, u, v
        real(dp) :: Li, Lj, Lk, dLi, dLj, dLk
        real(dp) :: Lu, Lv, dLu, dLv
        
        call timer_start('EVAL - Basis Functions')

        call finite_1D_positions(FE)

        if (allocated(FE%basis_at_quad))       deallocate(FE%basis_at_quad)
        if (allocated(FE%dbasis_dxi))  deallocate(FE%dbasis_dxi)
        if (allocated(FE%dbasis_deta)) deallocate(FE%dbasis_deta)
        if (allocated(FE%dbasis_dzeta)) deallocate(FE%dbasis_dzeta)
        if (allocated(FE%basis_at_bound_quad))     deallocate(FE%basis_at_bound_quad)
        if (allocated(FE%dbasis_bound_dxi))    deallocate(FE%dbasis_bound_dxi)
        if (allocated(FE%dbasis_bound_deta))   deallocate(FE%dbasis_bound_deta)
        if (allocated(FE%face_node_map)) deallocate(FE%face_node_map)

        allocate(FE%basis_at_quad(Quad3D%n_points, FE%n_basis))
        allocate(FE%dbasis_dxi(Quad3D%n_points, FE%n_basis))
        allocate(FE%dbasis_deta(Quad3D%n_points, FE%n_basis))
        allocate(FE%dbasis_dzeta(Quad3D%n_points, FE%n_basis))
        
        allocate(FE%basis_at_bound_quad(Quad2D%n_points, FE%n_nodes_per_face))
        allocate(FE%dbasis_bound_dxi(Quad2D%n_points, FE%n_nodes_per_face))
        allocate(FE%dbasis_bound_deta(Quad2D%n_points, FE%n_nodes_per_face))
        
        allocate(FE%face_node_map(FE%n_nodes_per_face, 6))

        ! 3D Hexahedral Face Map
        ! Faces: 1=Bottom(z=-1), 2=Top(z=+1), 3=Front(y=-1), 4=Right(x=+1), 5=Back(y=+1), 6=Left(x=-1)
        ! Basis index: n = i + (j-1)*(p+1) + (k-1)*(p+1)^2
        do n = 1, FE%n_nodes_per_face
            u = mod(n-1, FE%order+1) + 1
            v = (n-1)/(FE%order+1) + 1
            
            ! Face 1: Bottom (z = -1, k=1). Varies x(u), y(v)
            FE%face_node_map(n, 1) = u + (v-1)*(FE%order+1)

            ! Face 2: Top (z = +1, k=p+1). Varies x(u), y(v)
            FE%face_node_map(n, 2) = u + (v-1)*(FE%order+1) + (FE%order)*(FE%order+1)**2

            ! Face 3: Front (y = -1, j=1). Varies x(u), z(v)
            FE%face_node_map(n, 3) = u + (v-1)*(FE%order+1)**2

            ! Face 4: Right (x = +1, i=p+1). Varies y(u), z(v)
            FE%face_node_map(n, 4) = (FE%order+1) + (u-1)*(FE%order+1) + (v-1)*(FE%order+1)**2

            ! Face 5: Back (y = +1, j=p+1). Varies x(u), z(v)
            FE%face_node_map(n, 5) = u + (FE%order)*(FE%order+1) + (v-1)*(FE%order+1)**2

            ! Face 6: Left (x = -1, i=1). Varies y(u), z(v)
            FE%face_node_map(n, 6) = 1 + (u-1)*(FE%order+1) + (v-1)*(FE%order+1)**2
        end do

        do q = 1, Quad3D%n_points
            n = 1
            do k = 1, FE%order + 1
                Lk  = FE_basis_1D(FE, k, Quad3D%zeta(q))
                dLk = FE_basis_derivative_1D(FE, k, Quad3D%zeta(q))
                do j = 1, FE%order + 1
                    Lj  = FE_basis_1D(FE, j, Quad3D%eta(q))
                    dLj = FE_basis_derivative_1D(FE, j, Quad3D%eta(q))
                    do i = 1, FE%order + 1
                        Li  = FE_basis_1D(FE, i, Quad3D%xi(q))
                        dLi = FE_basis_derivative_1D(FE, i, Quad3D%xi(q))
                        
                        FE%basis_at_quad(q, n) = Li * Lj * Lk
                        FE%dbasis_dxi(q, n)    = dLi * Lj * Lk
                        FE%dbasis_deta(q, n)   = Li * dLj * Lk
                        FE%dbasis_dzeta(q, n)  = Li * Lj * dLk
                        n = n + 1
                    end do
                end do
            end do
        end do

        do q = 1, Quad2D%n_points
            n = 1
            do j = 1, FE%order + 1
                Lv = FE_basis_1D(FE, j, Quad2D%eta(q))
                dLv = FE_basis_derivative_1D(FE, j, Quad2D%eta(q))
                do i = 1, FE%order + 1
                    Lu = FE_basis_1D(FE, i, Quad2D%xi(q))
                    dLu = FE_basis_derivative_1D(FE, i, Quad2D%xi(q))
                    
                    FE%basis_at_bound_quad(q, n) = Lu * Lv
                    FE%dbasis_bound_dxi(q, n)    = dLu * Lv
                    FE%dbasis_bound_deta(q, n)   = Lu * dLv
                    n = n + 1
                end do
            end do
        end do

        call timer_stop('Basis Function EVAL')
    end subroutine InitialiseFiniteElements

    subroutine GetMapping(FE, q, elem_coords, dN_dx, dN_dy, dN_dz, detJ)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: q            
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: dN_dz(:)
        real(dp), intent(out)        :: detJ

        real(dp) :: J(3,3), invJ(3,3), dNdXiEta(3, FE%n_basis)
        integer  :: i ! For error message loop

        dNdXiEta(1, :) = FE%dbasis_dxi(q, :)
        dNdXiEta(2, :) = FE%dbasis_deta(q, :)
        dNdXiEta(3, :) = FE%dbasis_dzeta(q, :)
        J = matmul(dNdXiEta, elem_coords)

        detJ = J(1,1)*(J(2,2)*J(3,3) - J(2,3)*J(3,2)) - &
               J(1,2)*(J(2,1)*J(3,3) - J(2,3)*J(3,1)) + &
               J(1,3)*(J(2,1)*J(3,2) - J(2,2)*J(3,1))

        if (detJ <= 1.0d-10) then
            print *, " "
            print *, ">>> GEOMETRY ERROR: detJ =", detJ
            print *, ">>> Check node ordering for element coordinates below:"
            do i = 1, FE%n_basis
                write(*, '(A, I2, A, F10.5, A, F10.5, A, F10.5, A)') "Node ", i, ": (", &
                    elem_coords(i,1), ", ", elem_coords(i,2), ", ", elem_coords(i,3), ")"
            end do
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

        dN_dx = invJ(1,1) * FE%dbasis_dxi(q, :) + invJ(1,2) * FE%dbasis_deta(q, :) + invJ(1,3) * FE%dbasis_dzeta(q, :)
        dN_dy = invJ(2,1) * FE%dbasis_dxi(q, :) + invJ(2,2) * FE%dbasis_deta(q, :) + invJ(2,3) * FE%dbasis_dzeta(q, :)
        dN_dz = invJ(3,1) * FE%dbasis_dxi(q, :) + invJ(3,2) * FE%dbasis_deta(q, :) + invJ(3,3) * FE%dbasis_dzeta(q, :)

    end subroutine GetMapping

    subroutine GetArbitraryBasis(FE, xi, eta, zeta, N)
        type(t_finite), intent(in) :: FE
        real(dp), intent(in) :: xi, eta, zeta
        real(dp), intent(out) :: N(:)
        integer :: i, j, k, nn
        nn = 1
        do k = 1, FE%order + 1
            do j = 1, FE%order + 1
                do i = 1, FE%order + 1
                    N(nn) = FE_basis_1D(FE, i, xi) * FE_basis_1D(FE, j, eta) * FE_basis_1D(FE, k, zeta)
                    nn = nn + 1
                end do
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