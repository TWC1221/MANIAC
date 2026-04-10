module m_transport_precompute
    use m_types
    use m_quadrature
    use m_finite_elements
    use m_constants, only: dp_EPSILON, check_nan_scalar, check_nan_array, check_nan_matrix
    use m_material

    implicit none

    interface
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer, intent(in) :: m, n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
        end subroutine dgetrf

        subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            integer, intent(in) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgetrs
    end interface

contains

    subroutine InitialiseTransport(mesh, FE, Quad, QuadBound, sn_quad, materials, n_groups)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad, QuadBound
        type(t_sn_quadrature), intent(in) :: sn_quad
        type(t_material), intent(in) :: materials(:)
        integer, intent(in) :: n_groups

        ! 1. Geometric Integrals Mass and Stiffness
        call Precompute_Transport_Integrals(mesh, FE, Quad, QuadBound)
        
        ! 2. Angle Mapping 
        call Precompute_Reflective_Map(mesh, sn_quad)
        
        ! 3. LU Factorization
        call Precompute_All_Local_LU(mesh, FE, sn_quad, materials, n_groups)
        
    end subroutine InitialiseTransport

    subroutine Precompute_Transport_Integrals(mesh, FE, Quad, QuadBound)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad, QuadBound
        integer :: ee, q, ii, jj, f, i_f, j_f
        integer :: si, sj, n_cp_xi, n_cp_eta, p, n_basis, n_knots_xi, n_knots_eta
        real(dp) :: k_xi(size(mesh%knot_vectors_xi,2)), k_eta(size(mesh%knot_vectors_eta,2)), xi_val, eta_val, d_xi, d_eta
        real(dp) :: N_xi(FE%order+1), N_eta(FE%order+1), dN_xi(FE%order+1), dN_eta(FE%order+1), W, dW_dxi, dW_deta
        real(dp) :: nodes(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
        real(dp) :: coords(FE%n_nodes_per_face, 2), dx_dxi, dy_dxi
        real(dp) :: basis_local(FE%n_basis), weights_local(FE%n_basis), dN_dxi(FE%n_basis), dN_deta(FE%n_basis), w_face
        real(dp) :: J(2,2), invJ(2,2)

        n_knots_xi = size(mesh%knot_vectors_xi, 2)
        n_knots_eta = size(mesh%knot_vectors_eta, 2)
        allocate(mesh%elem_mass_matrix(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 mesh%elem_stiffness_x(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%elem_stiffness_y(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%face_mass_x(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems), &
                 mesh%face_mass_y(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems), &
                 mesh%basis_integrals_vol(FE%n_basis, mesh%n_elems))
        
        mesh%elem_mass_matrix = 0.0_dp; mesh%elem_stiffness_x = 0.0_dp; mesh%elem_stiffness_y = 0.0_dp
        mesh%face_mass_x = 0.0_dp; mesh%face_mass_y = 0.0_dp
        mesh%basis_integrals_vol = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(nodes, q, dN_dx, dN_dy, detJ, dV, ii, jj, f, coords, dx_dxi, dy_dxi, i_f, j_f, basis_local, &
        !$OMP& weights_local, si, sj, n_cp_xi, n_cp_eta, p, k_xi, k_eta, xi_val, eta_val, d_xi, d_eta, N_xi, N_eta, dN_xi, dN_eta, &
        !$OMP& W, dW_dxi, dW_deta, dN_dxi, dN_deta, J, invJ, n_basis)
        do ee = 1, mesh%n_elems
            if (any(mesh%elems(ee, 1:FE%n_basis) <= 0)) then
                ! Fallback: If connectivity is missing, we cannot integrate this element.
                cycle
            end if

            call check_nan_array(real(mesh%elems(ee, 1:FE%n_basis), dp), "mesh%elems", "Precompute_Transport_Integrals, Elem "//int_to_str(ee))
            nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), 1:2)
            
            ! Safety: Default weights to 1.0 if not provided or all zero
            if (.not. allocated(mesh%weights)) then
                weights_local = 1.0_dp
            else if (sum(abs(mesh%weights(mesh%elems(ee, 1:FE%n_basis)))) < 1.0d-12) then
                weights_local = 1.0_dp
            else
                weights_local = mesh%weights(mesh%elems(ee, 1:FE%n_basis))
            end if

            call check_nan_matrix(nodes, "nodes", "Precompute_Transport_Integrals, Elem "//int_to_str(ee))
            call check_nan_array(weights_local, "weights_local", "Precompute_Transport_Integrals, Elem "//int_to_str(ee))
            p = FE%order; n_cp_xi = p + 1; n_cp_eta = p + 1; n_basis = FE%n_basis
            k_xi = mesh%knot_vectors_xi(ee, :); k_eta = mesh%knot_vectors_eta(ee, :)
            call check_nan_array(k_xi, "k_xi", "Precompute_Transport_Integrals, Elem "//int_to_str(ee))
            
            ! Loop over non-zero knot spans in both parametric directions
            do sj = 1, n_knots_eta - 1
                d_eta = k_eta(sj+1) - k_eta(sj)
                if (d_eta < 1.0d-14) cycle
                do si = 1, n_knots_xi - 1
                    d_xi = k_xi(si+1) - k_xi(si)
                    if (d_xi < 1.0d-14) cycle
                    if (d_xi < dp_EPSILON) cycle

                    do q = 1, Quad%n_points
                        ! Map reference Gauss points to span-specific intervals
                        xi_val  = 0.5_dp * ((1.0_dp - Quad%xi(q)) * k_xi(si) + (1.0_dp + Quad%xi(q)) * k_xi(si+1))
                        eta_val = 0.5_dp * ((1.0_dp - Quad%eta(q)) * k_eta(sj) + (1.0_dp + Quad%eta(q)) * k_eta(sj+1))
                        call check_nan_scalar(xi_val, "xi_val", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_scalar(eta_val, "eta_val", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))

                        do ii = 1, n_cp_xi
                            N_xi(ii) = FE_basis_1D(k_xi, p, ii, xi_val)
                            dN_xi(ii) = FE_basis_derivative_1D(k_xi, p, ii, xi_val)
                        end do
                        do jj = 1, n_cp_eta
                            N_eta(jj) = FE_basis_1D(k_eta, p, jj, eta_val)
                            dN_eta(jj) = FE_basis_derivative_1D(k_eta, p, jj, eta_val)
                        end do
                        call check_nan_array(N_xi, "N_xi", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_array(N_eta, "N_eta", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_array(dN_xi, "dN_xi", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_array(dN_eta, "dN_eta", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))

                        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp
                        do jj = 1, n_cp_eta; do ii = 1, n_cp_xi
                            f = (jj-1)*n_cp_xi + ii
                            call check_nan_scalar(N_xi(ii), "N_xi(ii)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                            call check_nan_scalar(N_eta(jj), "N_eta(jj)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                            basis_local(f) = N_xi(ii) * N_eta(jj) * weights_local(f)
                            W = W + basis_local(f)
                            dN_dxi(f)  = dN_xi(ii) * N_eta(jj) * weights_local(f)
                            dN_deta(f) = N_xi(ii) * dN_eta(jj) * weights_local(f)
                            dW_dxi  = dW_dxi  + dN_dxi(f)
                            dW_deta = dW_deta + dN_deta(f)
                        end do; end do
                        call check_nan_scalar(W, "W", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_scalar(dW_dxi, "dW_dxi", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_scalar(dW_deta, "dW_deta", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        
                        if (W < 1.0d-15) cycle

                        dN_dxi  = (dN_dxi * W - basis_local * dW_dxi) / W**2
                        dN_deta = (dN_deta * W - basis_local * dW_deta) / W**2
                        basis_local = basis_local / W
                        call check_nan_array(dN_dxi, "dN_dxi (rational)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_array(dN_deta, "dN_deta (rational)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_array(basis_local, "basis_local (rational)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))

                        J(1,1) = dot_product(dN_dxi,  nodes(:,1)); J(1,2) = dot_product(dN_dxi,  nodes(:,2))
                        J(2,1) = dot_product(dN_deta, nodes(:,1)); J(2,2) = dot_product(dN_deta, nodes(:,2))
                        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
                        
                        if (abs(detJ) < 1.0d-16) then
                             cycle
                        end if
                        
                        invJ(1,1) = J(2,2)/detJ; invJ(1,2) = -J(1,2)/detJ
                        invJ(2,1) = -J(2,1)/detJ; invJ(2,2) = J(1,1)/detJ
                        call check_nan_matrix(invJ, "invJ", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))

                        dN_dx = invJ(1,1) * dN_dxi + invJ(2,1) * dN_deta
                        dN_dy = invJ(1,2) * dN_dxi + invJ(2,2) * dN_deta
                        call check_nan_array(dN_dx, "dN_dx", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_array(dN_dy, "dN_dy", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        
                        dV = abs(detJ) * Quad%weights(q) * 0.25_dp * d_xi * d_eta
                        call check_nan_scalar(dV, "dV", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))

                        mesh%elem_mass_matrix(:,:,ee)  = mesh%elem_mass_matrix(:,:,ee)  + outer_product(basis_local, basis_local) * dV
                        mesh%elem_stiffness_x(:,:,ee)  = mesh%elem_stiffness_x(:,:,ee)  + outer_product(dN_dx, basis_local) * dV
                        mesh%elem_stiffness_y(:,:,ee)  = mesh%elem_stiffness_y(:,:,ee)  + outer_product(dN_dy, basis_local) * dV
                        mesh%basis_integrals_vol(:,ee) = mesh%basis_integrals_vol(:,ee) + basis_local * dV
                    end do
                end do
            end do
            do f = 1, mesh%n_faces_per_elem
                if (any(mesh%elems(ee, FE%face_node_map(:, f)) <= 0)) cycle
                coords = mesh%nodes(mesh%elems(ee, FE%face_node_map(:, f)), 1:2)
                call check_nan_matrix(coords, "face_coords", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", Face "//int_to_str(f))
                do q = 1, QuadBound%n_points
                    ! Rational basis for face (NURBS 1D)
                    w_face = dot_product(FE%basis_at_bound_quad(q,:), weights_local(FE%face_node_map(:, f)))
                    call check_nan_scalar(w_face, "w_face", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", Face "//int_to_str(f)//", QP "//int_to_str(q))
                    
                    dx_dxi = dot_product(FE%dbasis_at_bound_quad(q,:), coords(:,1)); dy_dxi = dot_product(FE%dbasis_at_bound_quad(q,:), coords(:,2))
                    call check_nan_scalar(dx_dxi, "dx_dxi (face)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", Face "//int_to_str(f)//", QP "//int_to_str(q))
                    call check_nan_scalar(dy_dxi, "dy_dxi (face)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", Face "//int_to_str(f)//", QP "//int_to_str(q))
                    ! n*dS = (dy/dxi, -dx/dxi) * W
                    do i_f = 1, FE%n_nodes_per_face
                        do j_f = 1, FE%n_nodes_per_face
                            ii = FE%face_node_map(i_f, f); jj = FE%face_node_map(j_f, f)
                            W = (FE%basis_at_bound_quad(q,i_f) * weights_local(ii) / w_face) * &
                                (FE%basis_at_bound_quad(q,j_f) * weights_local(jj) / w_face)
                            call check_nan_scalar(W, "W (face basis)", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", Face "//int_to_str(f)//", QP "//int_to_str(q))
                            mesh%face_mass_x(ii,jj,f,ee) = mesh%face_mass_x(ii,jj,f,ee) + W * dy_dxi * QuadBound%weights(q)
                            mesh%face_mass_y(ii,jj,f,ee) = mesh%face_mass_y(ii,jj,f,ee) - W * dx_dxi * QuadBound%weights(q)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine Precompute_Transport_Integrals

    subroutine Precompute_Reflective_Map(mesh, sn_quad)
        type(t_mesh), intent(inout) :: mesh
        type(t_sn_quadrature), intent(in) :: sn_quad
        integer :: ee, f, mm, m_iter
        real(dp) :: normal(3), dir(3), ref_dir(3), max_dot, dprod

        allocate(mesh%reflect_map(sn_quad%n_angles, mesh%n_faces_per_elem, mesh%n_elems))
        mesh%reflect_map = 0

        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                normal(1:2) = mesh%face_normals(:, f, ee)
                normal(3) = 0.0_dp
                do mm = 1, sn_quad%n_angles
                    dir = sn_quad%dirs(mm, :)
                    ref_dir = dir - 2.0_dp * dot_product(dir, normal) * normal
                    max_dot = -1.0_dp
                    do m_iter = 1, sn_quad%n_angles
                        dprod = dot_product(ref_dir, sn_quad%dirs(m_iter, :))
                        if (dprod > max_dot) then
                            max_dot = dprod
                            mesh%reflect_map(mm, f, ee) = m_iter
                        end if
                    end do
                end do
            end do
        end do
    end subroutine Precompute_Reflective_Map

    subroutine Precompute_All_Local_LU(mesh, FE, sn_quad, materials, n_groups)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: sn_quad
        type(t_material), intent(in) :: materials(:)
        integer, intent(in) :: n_groups
        integer :: ee, mm, g, f, info
        real(dp) :: A(FE%n_basis, FE%n_basis), dir(2), o_n

        allocate(mesh%local_lu(FE%n_basis, FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups), &
                 mesh%local_pivots(FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups))
        mesh%local_lu = 0.0_dp; mesh%local_pivots = 0

        !$OMP PARALLEL DO PRIVATE(mm, dir, ee, g, A, f, o_n, info) SCHEDULE(DYNAMIC)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:2)
            do ee = 1, mesh%n_elems
                if (mesh%material_ids(ee) <= 0 .or. sum(abs(mesh%elem_mass_matrix(:,:,ee))) < 1.0d-15) cycle

                do g = 1, n_groups
                    call check_nan_matrix(mesh%elem_mass_matrix(:,:,ee), "elem_mass_matrix", "Precompute_All_Local_LU, Elem "//int_to_str(ee)//", Group "//int_to_str(g))
                    A = materials(mesh%material_ids(ee))%SigmaT(g) * mesh%elem_mass_matrix(:,:,ee) - &
                        (dir(1)*mesh%elem_stiffness_x(:,:,ee) + dir(2)*mesh%elem_stiffness_y(:,:,ee))
                    do f = 1, mesh%n_faces_per_elem
                        o_n = dot_product(dir, mesh%face_normals(:, f, ee))
                        if (o_n > 0.0_dp) A = A + (dir(1)*mesh%face_mass_x(:,:,f,ee) + dir(2)*mesh%face_mass_y(:,:,f,ee))
                    end do
                    call check_nan_matrix(A, "Matrix A before LU", "Precompute_All_Local_LU, Elem "//int_to_str(ee)//", Angle "//int_to_str(mm)//", Group "//int_to_str(g))
                    call dgetrf(FE%n_basis, FE%n_basis, A, FE%n_basis, mesh%local_pivots(:,ee,mm,g), info)
                    
                    if (info /= 0) then
                        write(*,*) "FATAL: LU Factorization failed (Singular Matrix) in Elem:", ee, " Angle:", mm
                        stop
                    end if

                    mesh%local_lu(:,:,ee,mm,g) = A
                end do
            end do
        end do
    end subroutine Precompute_All_Local_LU

end module m_transport_precompute