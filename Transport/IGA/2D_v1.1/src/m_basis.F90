module m_basis
    use m_constants
    use m_types
    implicit none
    public
        
    contains

    subroutine InitialiseBasis(FE, mesh)
        type(t_finite), intent(inout) :: FE
        type(t_mesh), intent(in) :: mesh
        integer :: ii, p, q

        FE%order = mesh%order
        p = FE%order; q = FE%order
        FE%p_order = p; FE%q_order = q
        FE%n_basis = (p + 1) * (q + 1)
        FE%n_nodes_per_face = p + 1

        if (allocated(FE%face_node_map)) deallocate(FE%face_node_map)
        allocate(FE%face_node_map(FE%n_nodes_per_face, 4))

        ! Face Node Mapping (CCW Loop for Outward Normals)
        do ii = 1, p + 1 ! ii is the local index along the face (1 to p+1)
            ! Face 1: Bottom (eta = -1, xi from -1 to 1)
            FE%face_node_map(ii, 1) = ii
            ! Face 2: Right (xi = 1, eta from -1 to 1)
            FE%face_node_map(ii, 2) = (ii - 1) * (p + 1) + (p + 1)
            ! Face 3: Top (eta = 1, xi from 1 to -1)
            FE%face_node_map(ii, 3) = (q + 1) * (p + 1) - (ii - 1)
            ! Face 4: Left (xi = -1, eta from 1 to -1)
            FE%face_node_map(ii, 4) = (p + 1) * (q + 1 - ii) + 1
        end do
    end subroutine InitialiseBasis

    subroutine GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, elem_coords, dN_dx, dN_dy, detJ, R, R_mat, xi_custom, eta_custom, J_out)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: ee
        integer, intent(in)          :: q
        type(t_mesh), intent(in)     :: mesh
        type(t_quadrature), intent(in) :: Quad
        real(dp), intent(in)         :: u1, u2, v1, v2
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: detJ, R(:)
        real(dp), optional, intent(out) :: R_mat(:,:)
        real(dp), optional, intent(in)  :: xi_custom, eta_custom
        real(dp), optional, intent(out) :: J_out(2,2)

        real(dp) :: J(2,2), invJ(2,2), dRdXiEta(2, size(R)), detJ_param, detJ_raw, detJ_abs
        real(dp) :: dR_dxi(size(R)), dR_deta(size(R))
        real(dp) :: xi, eta, xi_ref, eta_ref

        if (present(xi_custom) .and. present(eta_custom)) then
            xi_ref = xi_custom
            eta_ref = eta_custom
        else
            xi_ref = Quad%Xi(q)
            eta_ref = Quad%Eta(q)
        end if

        ! Map from reference [-1,1] to patch parametric space [u1, u2]
        xi  = 0.5_dp * ((u2 - u1) * xi_ref + (u2 + u1))
        eta = 0.5_dp * ((v2 - v1) * eta_ref + (v2 + v1))
        detJ_param = 0.25_dp * (u2 - u1) * (v2 - v1)

        ! Evaluate NURBS basis functions and parametric derivatives
        call EvalNURBS2D(FE, ee, mesh, xi, eta, R, dR_dxi, dR_deta)

        ! Physical mapping Jacobian
        dRdXiEta(1, :) = dR_dxi
        dRdXiEta(2, :) = dR_deta
        J = matmul(dRdXiEta, elem_coords)

        if (present(J_out)) J_out = J

        detJ_raw = (J(1,1)*J(2,2) - J(1,2)*J(2,1))
        detJ     = detJ_raw * detJ_param ! This is the physical Jacobian determinant
        detJ_abs = abs(detJ)

        ! Detect inverted or degenerate elements
        if (detJ < 0.0_dp) then
            write(*,*) "CRITICAL: Inverted Jacobian at Elem:", ee, " q:", q, " detJ:", detJ
            ! For now, we'll proceed, but this is a serious mesh error.
        end if
        if (detJ_abs < dp_EPSILON) then
            write(*,*) "CRITICAL: Degenerate Jacobian at Elem:", ee, " q:", q, " detJ:", detJ
            detJ_abs = dp_EPSILON ! Prevent division by zero, but indicates a bad element
        end if

        ! Compute inverse of physical Jacobian for chain rule
        invJ(1,1) =  J(2,2) / detJ_raw ! Note: detJ_raw can be negative for inverted elements
        invJ(1,2) = -J(1,2) / detJ_raw
        invJ(2,1) = -J(2,1) / detJ_raw
        invJ(2,2) =  J(1,1) / detJ_raw

        ! Chain rule: dR/dx = dR/du * du/dx + dR/dv * dv/dx
        dN_dx = (invJ(1,1) * dR_dxi + invJ(1,2) * dR_deta)
        dN_dy = (invJ(2,1) * dR_dxi + invJ(2,2) * dR_deta)

        if (present(R_mat)) R_mat = spread(R, dim=2, ncopies=FE%n_basis) * spread(R, dim=1, ncopies=FE%n_basis)
    end subroutine GetMapping2D

    subroutine EvalNURBS2D(FE, ee, mesh, xi, eta, R, dR_dxi, dR_deta)
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: ee
        type(t_mesh), intent(in) :: mesh
        real(dp), intent(in) :: xi, eta 
        real(dp), intent(out) :: R(:), dR_dxi(:), dR_deta(:)
        integer :: p, q, span_xi, span_eta, ii, jj, idx, p_idx
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1)
        real(dp) :: W, dW_dxi, dW_deta, w_ij
        real(dp) :: invW, invW2

        p = FE%p_order; q = FE%q_order; p_idx = mesh%elem_patch_id(ee)
        R = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp

        ! Slice knot vectors to handle varying lengths per patch
        call FindSpan(mesh%patch_n_cp_xi(p_idx)-1, p, xi, &
                      mesh%patch_knot_vectors_xi(p_idx, 1:mesh%patch_n_knots_xi(p_idx)), span_xi)
        call FindSpan(mesh%patch_n_cp_eta(p_idx)-1, q, eta, &
                      mesh%patch_knot_vectors_eta(p_idx, 1:mesh%patch_n_knots_eta(p_idx)), span_eta)
        call DersBasisFuns(span_xi, xi, p, 1, &
                           mesh%patch_knot_vectors_xi(p_idx, 1:mesh%patch_n_knots_xi(p_idx)), dN_xi)
        call DersBasisFuns(span_eta, eta, q, 1, &
                           mesh%patch_knot_vectors_eta(p_idx, 1:mesh%patch_n_knots_eta(p_idx)), dN_eta)

        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp

        idx = 0
        do jj = 1, q + 1
            do ii = 1, p + 1
                idx = idx + 1
                w_ij = mesh%weights(mesh%elems(ee, idx))
                W = W + dN_xi(1, ii) * dN_eta(1, jj) * w_ij
                dW_dxi = dW_dxi + dN_xi(2, ii) * dN_eta(1, jj) * w_ij
                dW_deta = dW_deta + dN_xi(1, ii) * dN_eta(2, jj) * w_ij
            end do
        end do

        invW = 1.0_dp / W
        invW2 = invW * invW

        idx = 0
        do jj = 1, q + 1
            do ii = 1, p + 1
                idx = idx + 1
                w_ij = mesh%weights(mesh%elems(ee, idx))
                R(idx) = (dN_xi(1, ii) * dN_eta(1, jj) * w_ij) * invW
                dR_dxi(idx) = ((dN_xi(2, ii) * dN_eta(1, jj) * w_ij) * W - (dN_xi(1, ii) * dN_eta(1, jj) * w_ij) * dW_dxi) * invW2
                dR_deta(idx) = ((dN_xi(1, ii) * dN_eta(2, jj) * w_ij) * W - (dN_xi(1, ii) * dN_eta(1, jj) * w_ij) * dW_deta) * invW2
            end do
        end do
    end subroutine EvalNURBS2D

    subroutine EvalNURBS1D(p, n_cp, knots, weights, xi, R, dR_dxi)
        integer, intent(in) :: p, n_cp
        real(dp), intent(in) :: knots(:), weights(:)
        real(dp), intent(in) :: xi
        real(dp), intent(out) :: R(:), dR_dxi(:)
        integer :: span, ii
        real(dp) :: dN(2, p+1)
        real(dp) :: W, dW_dxi, w_local
        real(dp) :: invW, invW2

        R = 0.0_dp; dR_dxi = 0.0_dp

        call FindSpan(n_cp-1, p, xi, knots, span)
        call DersBasisFuns(span, xi, p, 1, knots, dN)

        W = 0.0_dp; dW_dxi = 0.0_dp
        do ii = 1, p + 1
            w_local = weights(ii)
            W = W + dN(1, ii) * w_local
            dW_dxi = dW_dxi + dN(2, ii) * w_local
        end do

        invW = 1.0_dp / W
        invW2 = invW * invW

        do ii = 1, p + 1
            w_local = weights(ii)
            R(ii) = (dN(1, ii) * w_local) * invW
            dR_dxi(ii) = ((dN(2, ii) * w_local) * W - (dN(1, ii) * w_local) * dW_dxi) * invW2
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

    subroutine DersBasisFuns(ii, u, p, n, UU, ders)
        integer, intent(in) :: ii, p, n
        real(dp), intent(in) :: u, UU(:)
        real(dp), intent(out) :: ders(n+1, p+1)
        real(dp) :: ndu(p+1, p+1), left(p+1), right(p+1), saved, temp, d, a(2, p+1)
        integer :: jj, r, k, s1, s2, rk, pk, j1, j2
        ndu(1,1) = 1.0_dp
        do jj = 1, p
            left(jj+1) = u - UU(ii+1-jj); right(jj+1) = UU(ii+jj) - u; saved = 0.0_dp
            do r = 1, jj
                ndu(jj+1, r) = right(r+1) + left(jj-r+2)
                temp = ndu(r, jj) / ndu(jj+1, r)
                ndu(r, jj+1) = saved + right(r+1) * temp; saved = left(jj-r+2) * temp
            end do
            ndu(jj+1, jj+1) = saved
        end do
        ders(1, :) = ndu(:, p+1)
        do r = 0, p
            s1 = 0; s2 = 1; a(1, 1) = 1.0_dp
            do k = 1, n
                d = 0.0_dp; rk = r - k; pk = p - k
                if (r >= k) then; a(s2+1, 1) = a(s1+1, 1) / ndu(pk+2, rk+1); d = a(s2+1, 1) * ndu(rk+1, pk+1); end if
                j1 = merge(1, -rk+1, rk >= -1); j2 = merge(k-1, p-r, r <= p-k)
                do jj = j1, j2; a(s2+1, jj+1) = (a(s1+1, jj+1) - a(s1+1, jj)) / ndu(pk+2, rk+jj+1); d = d + a(s2+1, jj+1) * ndu(rk+jj+1, pk+1); end do
                if (r <= pk) then; a(s2+1, k+1) = -a(s1+1, k) / ndu(pk+2, r+1); d = d + a(s2+1, k+1) * ndu(r+1, pk+1); end if
                ders(k+1, r+1) = d; jj = s1; s1 = s2; s2 = jj
            end do
        end do
        r = p
        do k = 1, n; ders(k+1, :) = ders(k+1, :) * r; r = r * (p - k); end do
    end subroutine DersBasisFuns

end module m_basis
