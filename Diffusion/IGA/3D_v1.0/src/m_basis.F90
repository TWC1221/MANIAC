module m_basis
    use m_constants
    use m_quadrature
    use m_types
    implicit none
    public
        
    contains

    subroutine InitialiseBasis(FE)
        type(t_finite), intent(inout) :: FE
        integer :: i, idx, Poly

        Poly = FE%order
        FE%p_order = Poly
        FE%q_order = Poly
        FE%r_order = Poly
        ! n_basis is patch-dependent and taken as max CP count in the global mesh

        if (allocated(FE%nodal_map)) deallocate(FE%nodal_map)
        allocate(FE%nodal_map(FE%n_basis))

        idx = 1
        do i = 1, FE%n_basis
            FE%nodal_map(i) = i
        end do

    end subroutine InitialiseBasis

    subroutine GetMapping2D(FE, ee, mesh, qq, Quad, u1, u2, v1, v2, elem_coords, dN_dx, dN_dy, detJ, R, R_mat)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: ee
        integer, intent(in)          :: qq
        type(t_mesh), intent(in)     :: mesh
        type(t_quadrature), intent(in) :: Quad
        real(dp), intent(in)         :: u1, u2, v1, v2
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:)
        real(dp), intent(out)        :: detJ, R(:)
        real(dp), optional, intent(out) :: R_mat(:,:)

        real(dp) :: J(2,3), dRdXiEta(2, size(R)), detJ_param
        real(dp) :: g1(3), g2(3), g11, g22, g12
        real(dp) :: dR_dxi(size(R)), dR_deta(size(R))
        real(dp) :: xi, eta
        integer  :: n_en

        ! Map from reference [-1,1] to knot span [u1, u2]
        xi = 0.5_dp * ((u2 - u1) * Quad%Xi(qq) + (u2 + u1))
        eta = 0.5_dp * ((v2 - v1) * Quad%Eta(qq) + (v2 + v1))
        detJ_param = 0.25_dp * (u2 - u1) * (v2 - v1)

        ! Evaluate NURBS basis functions and parametric derivatives
        call EvalNURBS2D(FE, ee, mesh, xi, eta, R, dR_dxi, dR_deta)

        n_en = size(R)
        dN_dx = 0.0_dp
        dN_dy = 0.0_dp

        ! Physical mapping Jacobian
        dRdXiEta(1, :) = dR_dxi
        dRdXiEta(2, :) = dR_deta
        J = matmul(dRdXiEta, elem_coords)

        ! For 3D boundaries, detJ is the magnitude of the normal vector (surface area element)
        g1 = J(1,:); g2 = J(2,:)
        g11 = dot_product(g1, g1); g22 = dot_product(g2, g2); g12 = dot_product(g1, g2)
        detJ = sqrt(max(0.0_dp, g11*g22 - g12**2)) * detJ_param

        if (present(R_mat)) R_mat = spread(R, dim=2, ncopies=n_en) * spread(R, dim=1, ncopies=n_en)
    end subroutine GetMapping2D

    subroutine GetMapping3D(FE, ee, mesh, qq, Quad, u1, u2, v1, v2, w1, w2, elem_coords, dN_dx, dN_dy, dN_dz, detJ, R, R_mat)
        type(t_finite), intent(in)   :: FE
        integer, intent(in)          :: ee
        integer, intent(in)          :: qq
        type(t_mesh), intent(in)     :: mesh
        type(t_quadrature), intent(in) :: Quad
        real(dp), intent(in)         :: u1, u2, v1, v2, w1, w2
        real(dp), intent(in)         :: elem_coords(:,:) 
        real(dp), intent(out)        :: dN_dx(:), dN_dy(:), dN_dz(:)
        real(dp), intent(out)        :: detJ, R(:)
        real(dp), optional, intent(out) :: R_mat(:,:)

        real(dp) :: J(3,3), invJ(3,3), dRdXiEtaZeta(3, size(R)), detJ_param, detJ_phys
        real(dp) :: dR_dxi(size(R)), dR_deta(size(R)), dR_dzeta(size(R))
        real(dp) :: xi, eta, zeta
        integer  :: n_en

        xi   = 0.5_dp * ((u2 - u1) * Quad%Xi(qq)  + (u2 + u1))
        eta  = 0.5_dp * ((v2 - v1) * Quad%Eta(qq) + (v2 + v1))
        zeta = 0.5_dp * ((w2 - w1) * Quad%Zeta(qq) + (w2 + w1))
        detJ_param = 0.125_dp * (u2 - u1) * (v2 - v1) * (w2 - w1)

        call EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R, dR_dxi, dR_deta, dR_dzeta)

        n_en = size(R)
        dN_dx = 0.0_dp; dN_dy = 0.0_dp; dN_dz = 0.0_dp

        dRdXiEtaZeta(1, :) = dR_dxi
        dRdXiEtaZeta(2, :) = dR_deta
        dRdXiEtaZeta(3, :) = dR_dzeta
        J = matmul(dRdXiEtaZeta, elem_coords)

        detJ_phys = (J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) - &
                     J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) + &
                     J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1)))

        invJ(1,1) = (J(2,2)*J(3,3)-J(2,3)*J(3,2))/detJ_phys
        invJ(1,2) = (J(1,3)*J(3,2)-J(1,2)*J(3,3))/detJ_phys
        invJ(1,3) = (J(1,2)*J(2,3)-J(1,3)*J(2,2))/detJ_phys
        invJ(2,1) = (J(2,3)*J(3,1)-J(2,1)*J(3,3))/detJ_phys
        invJ(2,2) = (J(1,1)*J(3,3)-J(1,3)*J(3,1))/detJ_phys
        invJ(2,3) = (J(1,3)*J(2,1)-J(1,1)*J(2,3))/detJ_phys
        invJ(3,1) = (J(2,1)*J(3,2)-J(2,2)*J(3,1))/detJ_phys
        invJ(3,2) = (J(1,2)*J(3,1)-J(1,1)*J(3,2))/detJ_phys
        invJ(3,3) = (J(1,1)*J(2,2)-J(1,2)*J(2,1))/detJ_phys

        ! Apply Chain Rule across the full basis
        dN_dx = invJ(1,1) * dR_dxi + invJ(1,2) * dR_deta + invJ(1,3) * dR_dzeta
        dN_dy = invJ(2,1) * dR_dxi + invJ(2,2) * dR_deta + invJ(2,3) * dR_dzeta
        dN_dz = invJ(3,1) * dR_dxi + invJ(3,2) * dR_deta + invJ(3,3) * dR_dzeta

        detJ = detJ_phys * detJ_param

        if (present(R_mat)) R_mat = spread(R, dim=2, ncopies=n_en) * spread(R, dim=1, ncopies=n_en)
    end subroutine GetMapping3D

    subroutine EvalNURBS2D(FE, ee, mesh, xi, eta, R, dR_dxi, dR_deta)
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: ee
        type(t_mesh), intent(in) :: mesh
        real(dp), intent(in) :: xi, eta 
        real(dp), intent(out) :: R(:), dR_dxi(:), dR_deta(:)
        integer :: pp, qq, span_xi, span_eta, i, j, node_id, Ni, Nj, loc_idx_global_cp, i_g, j_g
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1)
        real(dp) :: W, dW_dxi, dW_deta, w_ij
        real(dp) :: invW, invW2

        pp = FE%p_order; qq = FE%q_order
        Ni = mesh%n_cp_xi_edge(ee) ! Control points in Xi
        Nj = mesh%n_cp_eta_edge(ee) ! Control points in Eta

        ! Initialize outputs
        R = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp

        call FindSpan(Ni-1, pp, xi, mesh%edge_knots_xi(ee,:), span_xi)
        call FindSpan(Nj-1, qq, eta, mesh%edge_knots_eta(ee,:), span_eta)
        call DersBasisFuns(span_xi, xi, pp, 1, mesh%edge_knots_xi(ee,:), dN_xi)
        call DersBasisFuns(span_eta, eta, qq, 1, mesh%edge_knots_eta(ee,:), dN_eta)

        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp
        do j = 1, qq + 1
            j_g = span_eta - qq + j - 2 ! 0-indexed global basis function index for Eta
            do i = 1, pp + 1
                i_g = span_xi - pp + i - 2 ! 0-indexed global basis function index for Xi
                loc_idx_global_cp = j_g * Ni + i_g + 1
                node_id = mesh%edges(ee, loc_idx_global_cp)
                w_ij = mesh%weights(node_id)
                
                W = W + dN_xi(1, i) * dN_eta(1, j) * w_ij
                dW_dxi = dW_dxi + dN_xi(2, i) * dN_eta(1, j) * w_ij
                dW_deta = dW_deta + dN_xi(1, i) * dN_eta(2, j) * w_ij
            end do
        end do
        if (abs(W) < 1e-16_dp) W = 1.0_dp
        invW = 1.0_dp / W
        invW2 = invW * invW

        do j = 1, qq + 1
            j_g = span_eta - qq + j - 2
            do i = 1, pp + 1
                i_g = span_xi - pp + i - 2
                loc_idx_global_cp = j_g * Ni + i_g + 1
                node_id = mesh%edges(ee, loc_idx_global_cp)
                w_ij = merge(mesh%weights(node_id), 1.0_dp, node_id > 0)
                R(loc_idx_global_cp) = (dN_xi(1, i) * dN_eta(1, j) * w_ij) * invW
                dR_dxi(loc_idx_global_cp) = ( (dN_xi(2, i) * dN_eta(1, j) * w_ij) * W - (dN_xi(1, i) * dN_eta(1, j) * w_ij) * dW_dxi ) * invW2
                dR_deta(loc_idx_global_cp) = ( (dN_xi(1, i) * dN_eta(2, j) * w_ij) * W - (dN_xi(1, i) * dN_eta(1, j) * w_ij) * dW_deta ) * invW2
            end do
        end do
    end subroutine EvalNURBS2D

    subroutine EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R, dR_dxi, dR_deta, dR_dzeta)
        type(t_finite), intent(in) :: FE
        integer, intent(in)        :: ee
        type(t_mesh), intent(in)   :: mesh
        real(dp), intent(in)       :: xi, eta, zeta
        real(dp), intent(out)      :: R(:), dR_dxi(:), dR_deta(:), dR_dzeta(:)
        
        integer :: pp, qq, rr, span_xi, span_eta, span_zeta, i, j, k, node_id
        integer :: Ni, Nj, Nk, loc_idx_global_cp, j_g, i_g, k_g
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1), dN_zeta(2, FE%r_order+1)
        real(dp) :: W, dW_dxi, dW_deta, dW_dzeta, w_ijk, Ri_num
        real(dp) :: invW, invW2

        pp = FE%p_order; qq = FE%q_order; rr = FE%r_order
        Ni = mesh%n_cp_xi(ee)
        Nj = mesh%n_cp_eta(ee)
        Nk = mesh%n_cp_zeta(ee)

        ! Initialize outputs
        R = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp; dR_dzeta = 0.0_dp

        ! 1. Find the knot spans
        call FindSpan(Ni-1, pp, xi,   mesh%knot_vectors_xi(ee,:),   span_xi)
        call FindSpan(Nj-1, qq, eta,  mesh%knot_vectors_eta(ee,:),  span_eta)
        call FindSpan(Nk-1, rr, zeta, mesh%knot_vectors_zeta(ee,:), span_zeta)

        ! 2. Evaluate B-spline basis functions and their 1st derivatives
        call DersBasisFuns(span_xi,   xi,   pp, 1, mesh%knot_vectors_xi(ee,:),   dN_xi)
        call DersBasisFuns(span_eta,  eta,  qq, 1, mesh%knot_vectors_eta(ee,:),  dN_eta)
        call DersBasisFuns(span_zeta, zeta, rr, 1, mesh%knot_vectors_zeta(ee,:), dN_zeta)

        ! 3. First Pass: Compute the NURBS denominator (W) and its derivatives
        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp; dW_dzeta = 0.0_dp
        
        do j = 1, qq + 1
            j_g = span_eta - qq + j - 2
            do i = 1, pp + 1
                i_g = span_xi - pp + i - 2
                do k = 1, rr + 1
                    k_g = span_zeta - rr + k - 2
                    ! Mapping: Eta (slowest) -> Xi (middle) -> Zeta (fastest)
                    loc_idx_global_cp = j_g * (Ni * Nk) + i_g * Nk + k_g + 1
                    node_id = mesh%elems(ee, loc_idx_global_cp)

                    w_ijk = merge(mesh%weights(node_id), 1.0_dp, node_id > 0)

                    W        = W        + dN_xi(1, i) * dN_eta(1, j) * dN_zeta(1, k) * w_ijk
                    dW_dxi   = dW_dxi   + dN_xi(2, i) * dN_eta(1, j) * dN_zeta(1, k) * w_ijk
                    dW_deta  = dW_deta  + dN_xi(1, i) * dN_eta(2, j) * dN_zeta(1, k) * w_ijk
                    dW_dzeta = dW_dzeta + dN_xi(1, i) * dN_eta(1, j) * dN_zeta(2, k) * w_ijk
                end do
            end do
        end do

        ! 4. Second Pass: Compute NURBS basis functions and their spatial derivatives
        invW  = 1.0_dp / W
        invW2 = invW * invW
        do j = 1, qq + 1
            j_g = span_eta - qq + j - 2
            do i = 1, pp + 1
                i_g = span_xi - pp + i - 2
                do k = 1, rr + 1
                    k_g = span_zeta - rr + k - 2
                    loc_idx_global_cp = j_g * (Ni * Nk) + i_g * Nk + k_g + 1
                    node_id = mesh%elems(ee, loc_idx_global_cp)
                    w_ijk = merge(mesh%weights(node_id), 1.0_dp, node_id > 0)
                    
                    Ri_num = dN_xi(1, i) * dN_eta(1, j) * dN_zeta(1, k) * w_ijk
                    R(loc_idx_global_cp) = Ri_num * invW
                    dR_dxi(loc_idx_global_cp) = ((dN_xi(2, i) * dN_eta(1, j) * dN_zeta(1, k) * w_ijk) * W - Ri_num * dW_dxi) * invW2
                    dR_deta(loc_idx_global_cp) = ((dN_xi(1, i) * dN_eta(2, j) * dN_zeta(1, k) * w_ijk) * W - Ri_num * dW_deta) * invW2
                    dR_dzeta(loc_idx_global_cp) = ((dN_xi(1, i) * dN_eta(1, j) * dN_zeta(2, k) * w_ijk) * W - Ri_num * dW_dzeta) * invW2
                end do
            end do
        end do
    
    end subroutine EvalNURBS3D

    subroutine FindSpan(n, pp, u, UU, span)
        integer, intent(in) :: n, pp
        real(dp), intent(in) :: u, UU(:)
        integer, intent(out) :: span
        integer :: low, high, mid

        ! Safety check to prevent out-of-bounds access if Ni was read incorrectly
        if (n + 2 > size(UU)) then; span = pp + 1; return; end if

        if (u >= UU(n+2)) then
            span = n + 1
            return
        end if
        if (u <= UU(pp+1)) then
            span = pp + 1
            return
        end if

        low = pp + 1; high = n + 2; mid = (low + high) / 2
        do while (u < UU(mid) .or. u >= UU(mid+1))
            if (u < UU(mid)) then; high = mid; else; low = mid; end if
            mid = (low + high) / 2
        end do
        span = mid
    end subroutine FindSpan

    subroutine DersBasisFuns(i, u, pp, n, UU, ders)
        integer, intent(in) :: i, pp, n
        real(dp), intent(in) :: u, UU(:)
        real(dp), intent(out) :: ders(n+1, pp+1)
        real(dp) :: ndu(pp+1, pp+1), left(pp+1), right(pp+1), saved, temp, d, a(2, pp+1)
        integer :: j, rr, k, s1, s2, rk, pk, j1, j2
        ndu(1,1) = 1.0_dp
        do j = 1, pp
            left(j+1) = u - UU(i+1-j); right(j+1) = UU(i+j) - u; saved = 0.0_dp
            do rr = 1, j
                ndu(j+1, rr) = right(rr+1) + left(j-rr+2)
                temp = ndu(rr, j) / ndu(j+1, rr)
                ndu(rr, j+1) = saved + right(rr+1) * temp; saved = left(j-rr+2) * temp
            end do
            ndu(j+1, j+1) = saved
        end do
        ders(1, :) = ndu(:, pp+1)
        do rr = 0, pp
            s1 = 0; s2 = 1; a(1, 1) = 1.0_dp
            do k = 1, n
                d = 0.0_dp; rk = rr - k; pk = pp - k
                if (rr >= k) then; a(s2+1, 1) = a(s1+1, 1) / ndu(pk+2, rk+1); d = a(s2+1, 1) * ndu(rk+1, pk+1); end if
                j1 = merge(1, -rk+1, rk >= -1); j2 = merge(k-1, pp-rr, rr <= pp-k)
                do j = j1, j2; a(s2+1, j+1) = (a(s1+1, j+1) - a(s1+1, j)) / ndu(pk+2, rk+j+1); d = d + a(s2+1, j+1) * ndu(rk+j+1, pk+1); end do
                if (rr <= pk) then; a(s2+1, k+1) = -a(s1+1, k) / ndu(pk+2, rr+1); d = d + a(s2+1, k+1) * ndu(rr+1, pk+1); end if
                ders(k+1, rr+1) = d; j = s1; s1 = s2; s2 = j
            end do
        end do
        rr = pp
        do k = 1, n; ders(k+1, :) = ders(k+1, :) * rr; rr = rr * (pp - k); end do
    end subroutine DersBasisFuns

end module m_basis