module m_basis
    use m_constants
    use m_types
    implicit none
    public
        
    contains

    subroutine InitialiseBasis(FE, mesh)
        type(t_finite), intent(inout) :: FE
        type(t_mesh), intent(in) :: mesh
        integer :: ii, jj, kk, p, q, r 

        FE%order = mesh%order
        p = FE%order; q = FE%order; r = FE%order
        FE%p_order = p; FE%q_order = q; FE%r_order = r
        FE%n_basis = (p + 1) * (q + 1) * (r + 1)
        FE%n_nodes_per_face = (p + 1) * (q + 1)

        allocate(FE%face_node_map(FE%n_nodes_per_face, 6))

        ! Face Node Mapping for Hexahedral Elements
        ! Local indexing: (k-1)*(p+1)*(q+1) + (j-1)*(p+1) + i
        
        ! Face 1: Bottom (zeta = -1, k = 1)
        do jj = 1, q + 1
            do ii = 1, p + 1
                FE%face_node_map((jj-1)*(p+1) + ii, 1) = (jj-1)*(p+1) + ii
            end do
        end do

        ! Face 2: Top (zeta = 1, k = r + 1)
        do jj = 1, q + 1
            do ii = 1, p + 1
                FE%face_node_map((jj-1)*(p+1) + ii, 2) = r*(p+1)*(q+1) + (jj-1)*(p+1) + ii
            end do
        end do

        ! Face 3: Front (eta = -1, j = 1)
        do kk = 1, r + 1
            do ii = 1, p + 1
                FE%face_node_map((kk-1)*(p+1) + ii, 3) = (kk-1)*(p+1)*(q+1) + ii
            end do
        end do

        ! Face 4: Back (eta = 1, j = q + 1)
        do kk = 1, r + 1
            do ii = 1, p + 1
                FE%face_node_map((kk-1)*(p+1) + ii, 4) = (kk-1)*(p+1)*(q+1) + q*(p+1) + ii
            end do
        end do

        ! Face 5: Left (xi = -1, i = 1)
        do kk = 1, r + 1
            do jj = 1, q + 1
                FE%face_node_map((kk-1)*(q+1) + jj, 5) = (kk-1)*(p+1)*(q+1) + (jj-1)*(p+1) + 1
            end do
        end do

        ! Face 6: Right (xi = 1, i = p + 1)
        do kk = 1, r + 1
            do jj = 1, q + 1
                FE%face_node_map((kk-1)*(q+1) + jj, 6) = (kk-1)*(p+1)*(q+1) + (jj-1)*(p+1) + p + 1
            end do
        end do

    end subroutine InitialiseBasis

    subroutine GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, elem_coords, dN_dx, dN_dy, detJ, R_basis, R_mat, xi_custom, eta_custom, J_out)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: ee
        integer, intent(in)          :: q
        type(t_mesh), intent(in)     :: mesh
        type(t_quadrature), intent(in) :: Quad
        real(dp), intent(in)         :: u1, u2, v1, v2
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: detJ, R_basis(:)
        real(dp), optional, intent(out) :: R_mat(:,:)
        real(dp), optional, intent(in)  :: xi_custom, eta_custom
        real(dp), optional, intent(out) :: J_out(2,2)

        real(dp) :: J(2,size(elem_coords,2)), invJ(2,2), dRdXiEta(2, size(R_basis)), detJ_param, detJ_raw, detJ_abs
        real(dp) :: dR_dxi(size(R_basis)), dR_deta(size(R_basis))
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
        call EvalNURBS2D(FE, ee, mesh, xi, eta, R_basis, dR_dxi, dR_deta)

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

        if (present(R_mat)) R_mat = spread(R_basis, dim=2, ncopies=FE%n_basis) * spread(R_basis, dim=1, ncopies=FE%n_basis)
    end subroutine GetMapping2D

    subroutine GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, elem_coords, dN_dx, dN_dy, dN_dz, detJ, R_basis, xi_custom, eta_custom, zeta_custom, J_out, R_out)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: ee, q
        type(t_mesh), intent(in)     :: mesh
        type(t_quadrature), intent(in) :: Quad
        real(dp), intent(in)         :: u1, u2, v1, v2, w1, w2
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:), dN_dz(:)
        real(dp), intent(out)        :: detJ, R_basis(:)
        real(dp), optional, intent(in)  :: xi_custom, eta_custom, zeta_custom
        real(dp), optional, intent(out) :: J_out(3,3)
        real(dp), optional, intent(out) :: R_out(3)

        real(dp) :: J(3,3), invJ(3,3), dRdXiEtaZeta(3, size(R_basis)), detJ_param, detJ_raw, detJ_abs
        real(dp) :: dR_dxi(size(R_basis)), dR_deta(size(R_basis)), dR_dzeta(size(R_basis))
        real(dp) :: xi, eta, zeta, xi_ref, eta_ref, zeta_ref

        if (present(xi_custom) .and. present(eta_custom) .and. present(zeta_custom)) then
            xi_ref = xi_custom; eta_ref = eta_custom; zeta_ref = zeta_custom
        else
            xi_ref = Quad%Xi(q); eta_ref = Quad%Eta(q); zeta_ref = Quad%Zeta(q)
        end if

        xi   = 0.5_dp * ((u2 - u1) * xi_ref + (u2 + u1))
        eta  = 0.5_dp * ((v2 - v1) * eta_ref + (v2 + v1))
        zeta = 0.5_dp * ((w2 - w1) * zeta_ref + (w2 + w1))
        detJ_param = 0.125_dp * (u2 - u1) * (v2 - v1) * (w2 - w1)

        call EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R_basis, dR_dxi, dR_deta, dR_dzeta)

        dRdXiEtaZeta(1, :) = dR_dxi
        dRdXiEtaZeta(2, :) = dR_deta
        dRdXiEtaZeta(3, :) = dR_dzeta
        J = matmul(dRdXiEtaZeta, elem_coords)

        if (present(J_out)) J_out = J

        if (present(R_out)) then
            R_out(1) = dot_product(R_basis, elem_coords(:, 1))
            R_out(2) = dot_product(R_basis, elem_coords(:, 2))
            R_out(3) = dot_product(R_basis, elem_coords(:, 3))
        end if

        detJ_raw = J(1,1)*(J(2,2)*J(3,3) - J(2,3)*J(3,2)) - &
                   J(1,2)*(J(2,1)*J(3,3) - J(2,3)*J(3,1)) + &
                   J(1,3)*(J(2,1)*J(3,2) - J(2,2)*J(3,1))
        
        detJ     = detJ_raw * detJ_param
        detJ_abs = abs(detJ)

        if (detJ < 0.0_dp) then
            write(*,*) "WARNING: Inverted 3D Jacobian at Elem:", ee, " q:", q, " detJ:", detJ
        end if
        if (detJ_abs < dp_EPSILON) then
            write(*,*) "CRITICAL: Degenerate 3D Jacobian at Elem:", ee
            detJ_abs = dp_EPSILON
        end if

        ! Adjugate matrix / detJ for inverse
        invJ(1,1) = (J(2,2)*J(3,3) - J(2,3)*J(3,2)) / detJ_raw
        invJ(1,2) = (J(1,3)*J(3,2) - J(1,2)*J(3,3)) / detJ_raw
        invJ(1,3) = (J(1,2)*J(2,3) - J(1,3)*J(2,2)) / detJ_raw
        invJ(2,1) = (J(2,3)*J(3,1) - J(2,1)*J(3,3)) / detJ_raw
        invJ(2,2) = (J(1,1)*J(3,3) - J(1,3)*J(3,1)) / detJ_raw
        invJ(2,3) = (J(1,3)*J(2,1) - J(1,1)*J(2,3)) / detJ_raw
        invJ(3,1) = (J(2,1)*J(3,2) - J(2,2)*J(3,1)) / detJ_raw
        invJ(3,2) = (J(1,2)*J(3,1) - J(1,1)*J(3,2)) / detJ_raw
        invJ(3,3) = (J(1,1)*J(2,2) - J(1,2)*J(2,1)) / detJ_raw

        dN_dx = (invJ(1,1) * dR_dxi + invJ(1,2) * dR_deta + invJ(1,3) * dR_dzeta)
        dN_dy = (invJ(2,1) * dR_dxi + invJ(2,2) * dR_deta + invJ(2,3) * dR_dzeta)
        dN_dz = (invJ(3,1) * dR_dxi + invJ(3,2) * dR_deta + invJ(3,3) * dR_dzeta)
    end subroutine GetMapping3D

    subroutine EvalNURBS2D(FE, ee, mesh, xi, eta, R_basis, dR_dxi, dR_deta)
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: ee
        type(t_mesh), intent(in) :: mesh
        real(dp), intent(in) :: xi, eta 
        real(dp), intent(out) :: R_basis(:), dR_dxi(:), dR_deta(:)
        integer :: p, q, span_xi, span_eta, ii, jj, idx, p_idx
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1)
        real(dp) :: W, dW_dxi, dW_deta, w_ij
        real(dp) :: invW, invW2

        p = FE%p_order; q = FE%q_order; p_idx = mesh%elem_patch_id(ee)
        R_basis = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp

        span_xi  = mesh%elem_span_indices(1, ee)
        span_eta = mesh%elem_span_indices(2, ee)
        
        call DersBasisFuns(span_xi, xi, p, 1, mesh%patches(p_idx)%knots_xi, dN_xi)
        call DersBasisFuns(span_eta, eta, q, 1, mesh%patches(p_idx)%knots_eta, dN_eta)

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
                R_basis(idx) = (dN_xi(1, ii) * dN_eta(1, jj) * w_ij) * invW
                dR_dxi(idx) = ((dN_xi(2, ii) * dN_eta(1, jj) * w_ij) * W - (dN_xi(1, ii) * dN_eta(1, jj) * w_ij) * dW_dxi) * invW2
                dR_deta(idx) = ((dN_xi(1, ii) * dN_eta(2, jj) * w_ij) * W - (dN_xi(1, ii) * dN_eta(1, jj) * w_ij) * dW_deta) * invW2
            end do
        end do
    end subroutine EvalNURBS2D

    subroutine EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R_basis, dR_dxi, dR_deta, dR_dzeta)
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: ee
        type(t_mesh), intent(in) :: mesh
        real(dp), intent(in) :: xi, eta, zeta
        real(dp), intent(out) :: R_basis(:), dR_dxi(:), dR_deta(:), dR_dzeta(:)
        integer :: p, q, r, span_xi, span_eta, span_zeta, ii, jj, kk, idx, p_idx
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1), dN_zeta(2, FE%r_order+1)
        real(dp) :: W, dW_dxi, dW_deta, dW_dzeta, w_ijk
        real(dp) :: invW, invW2

        p = FE%p_order; q = FE%q_order; r = FE%r_order
        p_idx = mesh%elem_patch_id(ee)
        R_basis = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp; dR_dzeta = 0.0_dp

        span_xi   = mesh%elem_span_indices(1, ee)
        span_eta  = mesh%elem_span_indices(2, ee)
        span_zeta = mesh%elem_span_indices(3, ee)

        call DersBasisFuns(span_xi, xi, p, 1, mesh%patches(p_idx)%knots_xi, dN_xi)
        call DersBasisFuns(span_eta, eta, q, 1, mesh%patches(p_idx)%knots_eta, dN_eta)
        call DersBasisFuns(span_zeta, zeta, r, 1, mesh%patches(p_idx)%knots_zeta, dN_zeta)

        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp; dW_dzeta = 0.0_dp

        idx = 0
        do kk = 1, r + 1
            do jj = 1, q + 1
                do ii = 1, p + 1
                    idx = idx + 1
                    w_ijk = mesh%weights(mesh%elems(ee, idx))
                    W = W + dN_xi(1, ii) * dN_eta(1, jj) * dN_zeta(1, kk) * w_ijk
                    dW_dxi = dW_dxi + dN_xi(2, ii) * dN_eta(1, jj) * dN_zeta(1, kk) * w_ijk
                    dW_deta = dW_deta + dN_xi(1, ii) * dN_eta(2, jj) * dN_zeta(1, kk) * w_ijk
                    dW_dzeta = dW_dzeta + dN_xi(1, ii) * dN_eta(1, jj) * dN_zeta(2, kk) * w_ijk
                end do
            end do
        end do

        invW = 1.0_dp / W
        invW2 = invW * invW

        idx = 0
        do kk = 1, r + 1
            do jj = 1, q + 1
                do ii = 1, p + 1
                    idx = idx + 1
                    w_ijk = mesh%weights(mesh%elems(ee, idx))
                    R_basis(idx) = (dN_xi(1, ii) * dN_eta(1, jj) * dN_zeta(1, kk) * w_ijk) * invW
                    dR_dxi(idx) = ((dN_xi(2, ii) * dN_eta(1, jj) * dN_zeta(1, kk) * w_ijk) * W - &
                                   (dN_xi(1, ii) * dN_eta(1, jj) * dN_zeta(1, kk) * w_ijk) * dW_dxi) * invW2
                    dR_deta(idx) = ((dN_xi(1, ii) * dN_eta(2, jj) * dN_zeta(1, kk) * w_ijk) * W - &
                                    (dN_xi(1, ii) * dN_eta(1, jj) * dN_zeta(1, kk) * w_ijk) * dW_deta) * invW2
                    dR_dzeta(idx) = ((dN_xi(1, ii) * dN_eta(1, jj) * dN_zeta(2, kk) * w_ijk) * W - &
                                     (dN_xi(1, ii) * dN_eta(1, jj) * dN_zeta(1, kk) * w_ijk) * dW_dzeta) * invW2
                end do
            end do
        end do
    end subroutine EvalNURBS3D

    subroutine FindSpan(n, p, u, UU, span)
        integer, intent(in) :: n, p
        real(dp), intent(in) :: u, UU(:)
        integer, intent(out) :: span
        integer :: low, high, mid
        if (u >= UU(n+2) - dp_EPSILON) then; span = n + 1; return; end if
        if (u <= UU(p+1)) then; span = p + 1; return; end if
        low = p + 1; high = n + 1
        do while (u < UU(mid) .or. u >= UU(mid+1))
            if (u < UU(mid)) then
                high = mid
            else
                low = mid
            end if
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
