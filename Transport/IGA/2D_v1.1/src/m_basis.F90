module m_basis
    use m_constants
    use m_quadrature
    use m_types
    implicit none
    public
        
    contains

    subroutine InitialiseBasis(FE, mesh)
        type(t_finite), intent(inout) :: FE
        type(t_mesh), intent(in) :: mesh
        integer :: i, p, ee, nx, ny

        FE%p_order = mesh%order
        FE%q_order = mesh%order
        p = mesh%order
        FE%n_nodes_per_face = max(maxval(mesh%n_cp_xi), maxval(mesh%n_cp_eta))
        FE%n_basis = mesh%nloc ! Synchronize basis with max control points per patch

        allocate(FE%face_node_map(FE%n_nodes_per_face, 4, mesh%n_elems))
        FE%face_node_map = 0

        do ee = 1, mesh%n_elems
            nx = mesh%n_cp_xi(ee)
            ny = mesh%n_cp_eta(ee)

            ! Face 1: Bottom (v_min) - Nodes 1 to nx
            do i = 1, nx; FE%face_node_map(i, 1, ee) = i; end do
            ! Face 2: Right (u_max) - Nodes nx, 2nx, ..., ny*nx
            do i = 1, ny; FE%face_node_map(i, 2, ee) = i * nx; end do
            ! Face 3: Top (v_max) - Nodes (ny-1)*nx + 1 to ny*nx
            do i = 1, nx; FE%face_node_map(i, 3, ee) = (ny - 1) * nx + i; end do
            ! Face 4: Left (u_min) - Nodes 1, nx+1, ..., (ny-1)*nx+1
            do i = 1, ny; FE%face_node_map(i, 4, ee) = (i - 1) * nx + 1; end do
        end do
    end subroutine InitialiseBasis

    subroutine GetMapping(FE, ee, mesh, q, Quad, u1, u2, v1, v2, nodes, dN_dx, dN_dy, detJ, R, R_mat)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: ee
        integer, intent(in)          :: q
        type(t_mesh), intent(in)     :: mesh
        type(t_quadrature), intent(in) :: Quad
        real(dp), intent(in)         :: u1, u2, v1, v2
        real(dp), intent(in)         :: nodes(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: detJ, R(:)
        real(dp), optional, intent(out) :: R_mat(:,:)

        real(dp) :: J(2,2), invJ(2,2), dRdXiEta(2, size(R)), detJ_param, detJ_phys
        real(dp) :: dR_dxi(size(R)), dR_deta(size(R))
        real(dp) :: xi, eta
        integer :: i

        ! Map from reference [-1,1] to knot span [u1, u2] x [v1, v2]
        xi = 0.5_dp * ((u2 - u1) * Quad%xi(q) + (u2 + u1))
        eta = 0.5_dp * ((v2 - v1) * Quad%eta(q) + (v2 + v1))
        detJ_param = 0.25_dp * (u2 - u1) * (v2 - v1)

        ! Evaluate NURBS basis functions and parametric derivatives
        call EvalNURBS2D(FE, ee, mesh, xi, eta, R, dR_dxi, dR_deta)

        ! Physical mapping Jacobian
        dRdXiEta(1, :) = dR_dxi
        dRdXiEta(2, :) = dR_deta
        J = matmul(dRdXiEta(:, 1:mesh%n_cp_xi(ee)*mesh%n_cp_eta(ee)), nodes(1:mesh%n_cp_xi(ee)*mesh%n_cp_eta(ee), :))

        detJ_phys = (J(1,1)*J(2,2) - J(1,2)*J(2,1))
        if (abs(detJ_phys) < 1e-15_dp) detJ_phys = sign(1e-15_dp, detJ_phys)

        detJ = detJ_phys * detJ_param

        invJ(1,1) =  J(2,2) / detJ_phys
        invJ(1,2) = -J(1,2) / detJ_phys
        invJ(2,1) = -J(2,1) / detJ_phys
        invJ(2,2) =  J(1,1) / detJ_phys

        ! Apply inverse Jacobian transpose for gradient mapping
        do i = 1, size(R)
            dN_dx(i) = invJ(1,1) * dR_dxi(i) + invJ(1,2) * dR_deta(i)
            dN_dy(i) = invJ(2,1) * dR_dxi(i) + invJ(2,2) * dR_deta(i)
        end do

        if (present(R_mat)) R_mat = spread(R, dim=2, ncopies=FE%n_basis) * spread(R, dim=1, ncopies=FE%n_basis)
    end subroutine GetMapping

    subroutine EvalNURBS2D(FE, ee, mesh, xi, eta, R, dR_dxi, dR_deta)
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: ee
        type(t_mesh), intent(in) :: mesh
        real(dp), intent(in) :: xi, eta 
        real(dp), intent(out) :: R(:), dR_dxi(:), dR_deta(:)
        integer :: p, q, span_xi, span_eta, i, j, idx
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1)
        real(dp) :: W, dW_dxi, dW_deta, w_ij
        real(dp) :: invW, invW2

        p = FE%p_order; q = FE%q_order
        R = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp

        call FindSpan(mesh%n_cp_xi(ee)-1, p, xi, mesh%knot_vectors_xi(ee, :), span_xi)
        call FindSpan(mesh%n_cp_eta(ee)-1, q, eta, mesh%knot_vectors_eta(ee, :), span_eta)
        call DersBasisFuns(span_xi, xi, p, 1, mesh%knot_vectors_xi(ee, :), dN_xi)
        call DersBasisFuns(span_eta, eta, q, 1, mesh%knot_vectors_eta(ee, :), dN_eta)

        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp
        do j = 1, q + 1
            do i = 1, p + 1
                idx = (span_eta - q + j - 2) * mesh%n_cp_xi(ee) + (span_xi - p + i - 1)

                w_ij = mesh%weights(mesh%elems(ee, idx))
                W = W + dN_xi(1, i) * dN_eta(1, j) * w_ij
                dW_dxi = dW_dxi + dN_xi(2, i) * dN_eta(1, j) * w_ij
                dW_deta = dW_deta + dN_xi(1, i) * dN_eta(2, j) * w_ij
            end do
        end do

        invW = 1.0_dp / W
        invW2 = invW * invW

        do j = 1, q + 1
            do i = 1, p + 1
                idx = (span_eta - q + j - 2) * mesh%n_cp_xi(ee) + (span_xi - p + i - 1)
                
                w_ij = mesh%weights(mesh%elems(ee, idx))
                R(idx) = (dN_xi(1, i) * dN_eta(1, j) * w_ij) * invW
                dR_dxi(idx) = ( (dN_xi(2, i) * dN_eta(1, j) * w_ij) * W - (dN_xi(1, i) * dN_eta(1, j) * w_ij) * dW_dxi ) * invW2
                dR_deta(idx) = ( (dN_xi(1, i) * dN_eta(2, j) * w_ij) * W - (dN_xi(1, i) * dN_eta(1, j) * w_ij) * dW_deta ) * invW2
            end do
        end do
    end subroutine EvalNURBS2D

    subroutine EvalNURBS1D(FE, ee, mesh, xi, R, dR_dxi)
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: ee ! edge index
        type(t_mesh), intent(in) :: mesh
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: R(:), dR_dxi(:)
        integer :: p, span, i, idx
        real(dp) :: dN(2, FE%p_order+1)
        real(dp) :: W, dW_dxi, w_i
        real(dp) :: invW, invW2

        p = FE%p_order
        R = 0.0_dp; dR_dxi = 0.0_dp

        call FindSpan(mesh%n_cp_edge(ee)-1, p, xi, mesh%edge_knots(ee, :), span)
        call DersBasisFuns(span, xi, p, 1, mesh%edge_knots(ee, :), dN)

        W = 0.0_dp; dW_dxi = 0.0_dp
        do i = 1, p + 1
            idx = span - p + i - 1
            w_i = mesh%weights(mesh%edges(ee, idx))
            W = W + dN(1, i) * w_i
            dW_dxi = dW_dxi + dN(2, i) * w_i
        end do

        invW = 1.0_dp / W
        invW2 = invW * invW

        do i = 1, p + 1
            idx = span - p + i - 1
            w_i = mesh%weights(mesh%edges(ee, idx))
            R(idx) = (dN(1, i) * w_i) * invW
            dR_dxi(idx) = ( (dN(2, i) * w_i) * W - (dN(1, i) * w_i) * dW_dxi ) * invW2
        end do
    end subroutine EvalNURBS1D

    subroutine FindSpan(n, p, u, UU, span)
        integer, intent(in) :: n, p
        real(dp), intent(in) :: u, UU(:)
        integer, intent(out) :: span
        integer :: low, high, mid
        if (u >= UU(n+2)) then; span = n + 1; return; end if
        if (u <= UU(p+1)) then; span = p + 1; return; end if
        low = p + 1; high = n + 1; mid = (low + high) / 2
        do while (u < UU(mid) .or. u >= UU(mid+1))
            if (u < UU(mid)) then; high = mid; else; low = mid; end if
            mid = (low + high) / 2
        end do
        span = mid
    end subroutine FindSpan

    subroutine DersBasisFuns(i, u, p, n, UU, ders)
        integer, intent(in) :: i, p, n
        real(dp), intent(in) :: u, UU(:)
        real(dp), intent(out) :: ders(n+1, p+1)
        real(dp) :: ndu(p+1, p+1), left(p+1), right(p+1), saved, temp, d, a(2, p+1)
        integer :: j, r, k, s1, s2, rk, pk, j1, j2
        ndu(1,1) = 1.0_dp
        do j = 1, p
            left(j+1) = u - UU(i+1-j); right(j+1) = UU(i+j) - u; saved = 0.0_dp
            do r = 1, j
                ndu(j+1, r) = right(r+1) + left(j-r+2)
                temp = ndu(r, j) / ndu(j+1, r)
                ndu(r, j+1) = saved + right(r+1) * temp; saved = left(j-r+2) * temp
            end do
            ndu(j+1, j+1) = saved
        end do
        ders(1, :) = ndu(:, p+1)
        do r = 0, p
            s1 = 0; s2 = 1; a(1, 1) = 1.0_dp
            do k = 1, n
                d = 0.0_dp; rk = r - k; pk = p - k
                if (r >= k) then; a(s2+1, 1) = a(s1+1, 1) / ndu(pk+2, rk+1); d = a(s2+1, 1) * ndu(rk+1, pk+1); end if
                j1 = merge(1, -rk+1, rk >= -1); j2 = merge(k-1, p-r, r <= p-k)
                do j = j1, j2; a(s2+1, j+1) = (a(s1+1, j+1) - a(s1+1, j)) / ndu(pk+2, rk+j+1); d = d + a(s2+1, j+1) * ndu(rk+j+1, pk+1); end do
                if (r <= pk) then; a(s2+1, k+1) = -a(s1+1, k) / ndu(pk+2, r+1); d = d + a(s2+1, k+1) * ndu(r+1, pk+1); end if
                ders(k+1, r+1) = d; j = s1; s1 = s2; s2 = j
            end do
        end do
        r = p
        do k = 1, n; ders(k+1, :) = ders(k+1, :) * r; r = r * (p - k); end do
    end subroutine DersBasisFuns

end module m_basis