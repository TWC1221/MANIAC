module m_transport_precompute
    use m_types
    use m_quadrature
    use m_basis
    use m_constants, only: dp_EPSILON, check_nan_scalar, check_nan_array, check_nan_matrix
    use m_material
    use m_sweep_order, only: get_face_cp_indices

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
        integer :: ee, q, ii, jj, f, i_f, j_f, ii_idx, jj_idx, n_patch_basis, i
        integer :: si, sj, p, n_basis, n_knots_xi, n_knots_eta
        real(dp) :: k_xi(size(mesh%knot_vectors_xi,2)), k_eta(size(mesh%knot_vectors_eta,2)), xi_val, eta_val, d_xi, d_eta, d_s
        real(dp) :: nodes(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
        real(dp) :: basis_local(FE%n_basis), weights_local(FE%n_basis), dN_dxi(FE%n_basis), dN_deta(FE%n_basis)
        real(dp) :: J(2,2), invJ(2,2), dx_ds, dy_ds, dW, fixed_val, val, fac
        integer :: idx_face(FE%n_basis), n_f_cp, s, range_start, range_end, count
        logical :: is_xi_fixed

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

        count = 0
        !$OMP PARALLEL DO PRIVATE(nodes, q, dN_dx, dN_dy, detJ, dV, ii, jj, f, i_f, j_f, basis_local, &
        !$OMP& weights_local, si, sj, p, k_xi, k_eta, xi_val, eta_val, d_xi, d_eta, d_s, i, n_patch_basis, &
        !$OMP& dN_dxi, dN_deta, J, invJ, n_basis, jj_idx, ii_idx, is_xi_fixed, fixed_val, range_start, range_end, idx_face, n_f_cp, fac, s, val, dx_ds, dy_ds, dW)
        do ee = 1, mesh%n_elems
            n_patch_basis = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee)
            if (n_patch_basis <= 0) cycle

            nodes = 0.0_dp; weights_local = 1.0_dp
            do i = 1, n_patch_basis
                if (mesh%elems(ee, i) > 0) then
                    nodes(i, :) = mesh%nodes(mesh%elems(ee, i), 1:2)
                    ! Ensure weights are never zero to avoid singular basis functions
                    if (mesh%weights(mesh%elems(ee, i)) < 1d-12) mesh%weights(mesh%elems(ee, i)) = 1.0_dp
                end if
            end do
            
            if (allocated(mesh%weights)) then
                if (sum(abs(weights_local(1:n_patch_basis))) < 1.0d-12) weights_local = 1.0_dp
            end if

            call check_nan_matrix(nodes, "nodes", "Precompute_Transport_Integrals, Elem "//int_to_str(ee))
            call check_nan_array(weights_local, "weights_local", "Precompute_Transport_Integrals, Elem "//int_to_str(ee))
            p = FE%order
            n_patch_basis = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee)
            k_xi = mesh%knot_vectors_xi(ee, :); k_eta = mesh%knot_vectors_eta(ee, :)

            ! Loop over non-zero knot spans in both parametric directions
            do sj = p + 1, mesh%n_cp_eta(ee)
                d_eta = k_eta(sj+1) - k_eta(sj)
                if (d_eta < 1.0d-14) cycle
                do si = p + 1, mesh%n_cp_xi(ee)
                    d_xi = k_xi(si+1) - k_xi(si)
                    if (d_xi < 1.0d-14) cycle

                    do q = 1, Quad%n_points
                        ! Map reference Gauss points to span-specific intervals
                        xi_val  = 0.5_dp * ((1.0_dp - Quad%xi(q)) * k_xi(si) + (1.0_dp + Quad%xi(q)) * k_xi(si+1))
                        eta_val = 0.5_dp * ((1.0_dp - Quad%eta(q)) * k_eta(sj) + (1.0_dp + Quad%eta(q)) * k_eta(sj+1))

                        call EvalNURBS2D(FE, ee, mesh, xi_val, eta_val, basis_local, dN_dxi, dN_deta)

                        J(1,1) = dot_product(dN_dxi,  nodes(:,1)); J(1,2) = dot_product(dN_dxi,  nodes(:,2))
                        J(2,1) = dot_product(dN_deta, nodes(:,1)); J(2,2) = dot_product(dN_deta, nodes(:,2))
                        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
                        if (abs(detJ) < 1.0e-18_dp) cycle
                        
                        if (detJ /= detJ .or. abs(detJ) < 1.0e-15_dp) cycle
                        
                        invJ(1,1) = J(2,2)/detJ;  invJ(1,2) = -J(1,2)/detJ
                        invJ(2,1) = -J(2,1)/detJ; invJ(2,2) =  J(1,1)/detJ
                        call check_nan_matrix(invJ, "invJ", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))

                        dN_dx = invJ(1,1) * dN_dxi + invJ(2,1) * dN_deta
                        dN_dy = invJ(1,2) * dN_dxi + invJ(2,2) * dN_deta

                        call check_nan_array(dN_dx, "dN_dx", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        call check_nan_array(dN_dy, "dN_dy", "Precompute_Transport_Integrals, Elem "//int_to_str(ee)//", QP "//int_to_str(q))
                        
                        dV = abs(detJ) * Quad%weights(q) * 0.25_dp * d_xi * d_eta
                        do jj_idx = 1, n_patch_basis
                            do ii_idx = 1, n_patch_basis
                                mesh%elem_mass_matrix(ii_idx, jj_idx, ee) = mesh%elem_mass_matrix(ii_idx, jj_idx, ee) + basis_local(ii_idx) * basis_local(jj_idx) * dV
                                mesh%elem_stiffness_x(ii_idx, jj_idx, ee) = mesh%elem_stiffness_x(ii_idx, jj_idx, ee) + dN_dx(ii_idx) * basis_local(jj_idx) * dV
                                mesh%elem_stiffness_y(ii_idx, jj_idx, ee) = mesh%elem_stiffness_y(ii_idx, jj_idx, ee) + dN_dy(ii_idx) * basis_local(jj_idx) * dV
                            end do
                            mesh%basis_integrals_vol(jj_idx, ee) = mesh%basis_integrals_vol(jj_idx, ee) + basis_local(jj_idx) * dV
                        end do
                    end do
                end do
            end do
            do f = 1, mesh%n_faces_per_elem
                select case(f)
                case(1); is_xi_fixed = .true.;  fixed_val = k_xi(1);  range_start = 1; range_end = mesh%n_knots_eta_patch(ee)
                case(2); is_xi_fixed = .true.;  fixed_val = k_xi(mesh%n_knots_xi_patch(ee)); range_start = 1; range_end = mesh%n_knots_eta_patch(ee)
                case(3); is_xi_fixed = .false.; fixed_val = k_eta(1); range_start = 1; range_end = mesh%n_knots_xi_patch(ee)
                case(4); is_xi_fixed = .false.; fixed_val = k_eta(mesh%n_knots_eta_patch(ee)); range_start = 1; range_end = mesh%n_knots_xi_patch(ee)
                end select
                
                call get_face_cp_indices(ee, f, mesh, idx_face, n_f_cp)
                fac = merge(-1.0_dp, 1.0_dp, f==1 .or. f==4)

                do s = range_start, range_end - 1
                    d_s = merge(k_eta(s+1)-k_eta(s), k_xi(s+1)-k_xi(s), is_xi_fixed)
                    if (d_s < 1.0d-14) cycle
                    do q = 1, QuadBound%n_points
                        val = 0.5_dp * ((1.0_dp - QuadBound%xi(q)) * merge(k_eta(s), k_xi(s), is_xi_fixed) + &
                                        (1.0_dp + QuadBound%xi(q)) * merge(k_eta(s+1), k_xi(s+1), is_xi_fixed))
                        if (is_xi_fixed) then
                            call EvalNURBS2D(FE, ee, mesh, fixed_val, val, basis_local, dN_dxi, dN_deta)
                            dx_ds = dot_product(dN_deta, nodes(:,1)); dy_ds = dot_product(dN_deta, nodes(:,2))
                        else
                            call EvalNURBS2D(FE, ee, mesh, val, fixed_val, basis_local, dN_dxi, dN_deta)
                            dx_ds = dot_product(dN_dxi, nodes(:,1)); dy_ds = dot_product(dN_dxi, nodes(:,2))
                        end if
                        if (dx_ds /= dx_ds .or. dy_ds /= dy_ds) cycle
                        
                        ! Use geometric normal components from the mesh object instead of manual 'fac'
                        dW = 0.5_dp * d_s * QuadBound%weights(q)

                        do jj_idx = 1, n_f_cp
                            jj = idx_face(jj_idx)
                            do ii = 1, n_patch_basis
                                mesh%face_mass_x(ii, jj, f, ee) = mesh%face_mass_x(ii, jj, f, ee) + basis_local(ii) * basis_local(jj) * mesh%face_normals(1, f, ee) * dW * sqrt(dx_ds**2 + dy_ds**2)
                                mesh%face_mass_y(ii, jj, f, ee) = mesh%face_mass_y(ii, jj, f, ee) + basis_local(ii) * basis_local(jj) * mesh%face_normals(2, f, ee) * dW * sqrt(dx_ds**2 + dy_ds**2)
                            end do
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
        integer :: ee, mm, g, f, info, n_patch_basis
        real(dp) :: A(FE%n_basis, FE%n_basis), dir(2), o_n, msum

        allocate(mesh%local_lu(FE%n_basis, FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups), &
                 mesh%local_pivots(FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups))
        mesh%local_lu = 0.0_dp; mesh%local_pivots = 0
        
        write(*,*) ">>> [DEBUG] Starting LU Precomputation. Checking for element stability."

        !$OMP PARALLEL DO PRIVATE(mm, dir, ee, g, A, f, o_n, info, n_patch_basis, msum) SCHEDULE(DYNAMIC)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:2)
            do ee = 1, mesh%n_elems
                n_patch_basis = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee)
                if (n_patch_basis <= 0) then
                    write(*,*) ">>> [DEBUG] Warning: Elem ", ee, " has no basis functions."
                    cycle
                end if

                do g = 1, n_groups
                    ! Ensure we have a valid material cross section, default to zero if ID is invalid
                    msum = 0.0_dp
                    if (mesh%material_ids(ee) > 0) msum = materials(mesh%material_ids(ee))%SigmaT(g)

                    A = msum * mesh%elem_mass_matrix(:,:,ee) - &
                        (dir(1)*mesh%elem_stiffness_x(:,:,ee) + dir(2)*mesh%elem_stiffness_y(:,:,ee))
                    
                    block
                        integer :: outflow_count
                        outflow_count = 0
                    do f = 1, mesh%n_faces_per_elem
                        o_n = dot_product(dir, mesh%face_normals(:, f, ee))
                        if (o_n > 1e-12_dp) then
                            A = A + (dir(1)*mesh%face_mass_x(:,:,f,ee) + dir(2)*mesh%face_mass_y(:,:,f,ee))
                            outflow_count = outflow_count + 1
                        end if
                    end do
                        if (outflow_count == 0) then
                            ! Critical diagnostic: Transport matrices are usually singular if there is no outflow
                            write(*,*) ">>> [DEBUG] WARNING: Elem ", ee, " Angle ", mm, " has ZERO outflow faces. Matrix will likely be singular."
                        end if
                    end block

                    call check_nan_matrix(A(1:n_patch_basis, 1:n_patch_basis), "Matrix A before LU", "Precompute_All_Local_LU, Elem "//int_to_str(ee)//", Angle "//int_to_str(mm)//", Group "//int_to_str(g))
                    call dgetrf(n_patch_basis, n_patch_basis, A, FE%n_basis, mesh%local_pivots(:,ee,mm,g), info)
                    
                    if (info /= 0) then
                        write(*,*) ">>> [DEBUG] FATAL: LU Factorization failed (Singular Matrix) in Elem:", ee, " Angle:", mm, " info:", info
                        write(*,*) ">>> [DEBUG] Diagonal of A for failed element:"
                        do f = 1, n_patch_basis
                            write(*,'(A,I4,A,ES12.4)') "   Diag(", f, "): ", A(f,f)
                        end do
                        write(*,*) ">>> [DEBUG] Elem Properties: SigmaT=", msum, " MassSum=", sum(mesh%elem_mass_matrix(:,:,ee))
                        write(*,*) ">>> [DEBUG] Direction vector: ", dir
                        stop
                    end if

                    mesh%local_lu(:,:,ee,mm,g) = A
                end do
            end do
        end do
    end subroutine Precompute_All_Local_LU

end module m_transport_precompute