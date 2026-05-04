module m_transport_precompute
    use m_types
    use m_quadrature
    use m_basis
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

    subroutine InitialiseTransport(mesh, FE, sn_quad, Quad2D, Quad1D, materials)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: sn_quad
        type(t_quadrature), intent(in) :: Quad2D, Quad1D
        type(t_material), intent(in) :: materials(:)
        
        ! 1. Precompute integrals (Matrices) - moved from main
        call Precompute_Transport_Integrals(mesh, FE, Quad2D, Quad1D)

        ! 2. Angle Mapping 
        call Precompute_Reflective_Map(mesh, sn_quad)
        
        ! 3. LU Factorization
        call Precompute_LU_Decomposition(mesh, FE, sn_quad, materials, mesh%n_groups)
        
        ! 4. Verification
        call Verify_Integrals(mesh, FE)

    end subroutine InitialiseTransport

    subroutine Precompute_Transport_Integrals(mesh, FE, Quad, QuadBound)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad, QuadBound
        integer :: ee, q, ii, jj, f
        real(dp) :: nodes(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
        real(dp) :: R(FE%n_basis), dx_dxi, dy_dxi, xi_f, eta_f, J(2,2)
        real(dp) :: u1, u2, v1, v2

        allocate(mesh%elem_mass_matrix(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 mesh%elem_stiffness_x(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%elem_stiffness_y(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%face_mass_x(FE%n_basis, FE%n_basis, 4, mesh%n_elems), &
                 mesh%face_mass_y(FE%n_basis, FE%n_basis, 4, mesh%n_elems), &
                 mesh%basis_integrals_vol(FE%n_basis, mesh%n_elems))
        
        mesh%elem_mass_matrix = 0.0_dp; mesh%elem_stiffness_x = 0.0_dp; mesh%elem_stiffness_y = 0.0_dp
        mesh%face_mass_x = 0.0_dp; mesh%face_mass_y = 0.0_dp
        mesh%basis_integrals_vol = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(ee, nodes, u1, u2, v1, v2, q, dN_dx, dN_dy, detJ, R, dV, f, xi_f, eta_f, J, dx_dxi, dy_dxi, ii, jj)
        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
            u1 = mesh%elem_u_min(ee); u2 = mesh%elem_u_max(ee)
            v1 = mesh%elem_v_min(ee); v2 = mesh%elem_v_max(ee)

            do q = 1, Quad%n_points
                call GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, nodes, dN_dx, dN_dy, detJ, R)
                dV = detJ * Quad%weights(q)
                mesh%elem_mass_matrix(:,:,ee)     = mesh%elem_mass_matrix(:,:,ee)     + spread(R, 2, FE%n_basis) * spread(R, 1, FE%n_basis) * dV
                mesh%elem_stiffness_x(:,:,ee)     = mesh%elem_stiffness_x(:,:,ee)     + spread(dN_dx, 2, FE%n_basis) * spread(R, 1, FE%n_basis) * dV
                mesh%elem_stiffness_y(:,:,ee)     = mesh%elem_stiffness_y(:,:,ee)     + spread(dN_dy, 2, FE%n_basis) * spread(R, 1, FE%n_basis) * dV
                mesh%basis_integrals_vol(:,ee)    = mesh%basis_integrals_vol(:,ee)    + R * dV
            end do
            
            do f = 1, 4 ! Number of faces per element
                do q = 1, QuadBound%n_points
                    ! Evaluate actual NURBS basis on the face coordinates
                    select case(f)
                        case(1); xi_f = QuadBound%xi(q); eta_f = -1.0_dp
                        case(2); xi_f = 1.0_dp;          eta_f = QuadBound%xi(q)
                        case(3); xi_f = -QuadBound%xi(q); eta_f = 1.0_dp
                        case(4); xi_f = -1.0_dp;         eta_f = -QuadBound%xi(q)
                    end select

                    call GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, nodes, dN_dx, dN_dy, detJ, R, &
                                     xi_custom=xi_f, eta_custom=eta_f, J_out=J)

                    ! Scale tangent components [dx/dxi_ref, dy/dxi_ref] for the parametric span.
                    ! Negation for faces 3 and 4 ensures a consistent CCW path for outward normals.
                    select case(f)
                        case(1); dx_dxi =  J(1,1)*0.5_dp*(u2-u1); dy_dxi =  J(1,2)*0.5_dp*(u2-u1)
                        case(2); dx_dxi =  J(2,1)*0.5_dp*(v2-v1); dy_dxi =  J(2,2)*0.5_dp*(v2-v1)
                        case(3); dx_dxi = -J(1,1)*0.5_dp*(u2-u1); dy_dxi = -J(1,2)*0.5_dp*(u2-u1)
                        case(4); dx_dxi = -J(2,1)*0.5_dp*(v2-v1); dy_dxi = -J(2,2)*0.5_dp*(v2-v1)
                    end select

                    do ii = 1, FE%n_basis
                        do jj = 1, FE%n_basis
                            mesh%face_mass_x(ii,jj,f,ee) = mesh%face_mass_x(ii,jj,f,ee) + R(ii) * R(jj) * dy_dxi * QuadBound%weights(q)
                            mesh%face_mass_y(ii,jj,f,ee) = mesh%face_mass_y(ii,jj,f,ee) - R(ii) * R(jj) * dx_dxi * QuadBound%weights(q)
                        end do
                    end do
                end do ! end q loop
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Precompute_Transport_Integrals

    subroutine Verify_Integrals(mesh, FE)
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        integer :: ee
        real(dp) :: vol, stiffness_err_x, stiffness_err_y, net_face_x, net_face_y, symmetry_err
        logical :: passed

        passed = .true.
        do ee = 1, mesh%n_elems
            ! 1. Volume Check
            vol = sum(mesh%basis_integrals_vol(:, ee))
            if (vol <= 0.0_dp) then
                write(*,*) "[ERROR] Element", ee, "has non-positive volume:", vol
                passed = .false.
            end if

            ! 2. Stiffness Sum Check (Weak Form: sum over test functions [dim=1] must be zero)
            stiffness_err_x = maxval(abs(sum(mesh%elem_stiffness_x(:, :, ee), dim=1)))
            stiffness_err_y = maxval(abs(sum(mesh%elem_stiffness_y(:, :, ee), dim=1)))

            ! 3. Divergence Theorem Check (Closed loop integral of normals = 0)
            net_face_x = sum(mesh%face_mass_x(:, :, :, ee))
            net_face_y = sum(mesh%face_mass_y(:, :, :, ee))

            ! 4. Symmetry Check for Mass Matrix
            symmetry_err = maxval(abs(mesh%elem_mass_matrix(:,:,ee) - transpose(mesh%elem_mass_matrix(:,:,ee))))

            if (max(stiffness_err_x, stiffness_err_y) > 1e-10_dp) then
                write(*,*) "[WARNING] Elem", ee, "Stiffness Consistency Error:", max(stiffness_err_x, stiffness_err_y)
            end if
            if (max(abs(net_face_x), abs(net_face_y)) > 1e-10_dp) then
                write(*,*) "[ERROR] Elem", ee, "Boundary Closure Error (Normal check):", net_face_x, net_face_y
                passed = .false.
            end if
            if (symmetry_err > 1e-12_dp) then
                write(*,*) "[WARNING] Elem", ee, "Mass Matrix Symmetry Error:", symmetry_err
            end if
        end do

        if (passed) write(*,*) ">>> Geometric Integrals Verified Successfully."
    end subroutine Verify_Integrals

    subroutine Precompute_Reflective_Map(mesh, sn_quad)
        type(t_mesh), intent(inout) :: mesh
        type(t_sn_quadrature), intent(in) :: sn_quad
        integer :: ee, f, mm, m_iter
        real(dp) :: normal(3), dir(3), ref_dir(3), max_dot, dprod

        allocate(mesh%reflect_map(sn_quad%n_angles, 4, mesh%n_elems))
        mesh%reflect_map = 0

        !$OMP PARALLEL DO PRIVATE(ee, f, normal, mm, dir, ref_dir, max_dot, m_iter, dprod)
        do ee = 1, mesh%n_elems
            do f = 1, 4 ! Number of faces per element
                normal(1:2) = mesh%face_normals(:, f, ee)
                normal(3) = 0.0_dp
                do mm = 1, sn_quad%n_angles
                    dir = sn_quad%dirs(mm, :)
                    ref_dir = dir - 2.0_dp * dot_product(dir, normal) * normal
                    
                    max_dot = -2.0_dp
                    do m_iter = 1, sn_quad%n_angles
                        ! STABILITY FIX: Reflections MUST stay on the same polar level (z-plane)
                        if (abs(ref_dir(3) - sn_quad%dirs(m_iter, 3)) > SMALL_NUMBER) cycle
                        
                        dprod = dot_product(ref_dir, sn_quad%dirs(m_iter, :))
                        if (dprod > max_dot) then
                            max_dot = dprod
                            mesh%reflect_map(mm, f, ee) = m_iter
                        end if
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Precompute_Reflective_Map

    subroutine Precompute_LU_Decomposition(mesh, FE, sn_quad, materials, n_groups)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: sn_quad
        type(t_material), intent(in) :: materials(:)
        integer, intent(in) :: n_groups
        integer :: ee, mm, g, f, info
        real(dp) :: A(FE%n_basis, FE%n_basis), dir(2), o_n
        real(dp) :: StiffnessOutflow(FE%n_basis, FE%n_basis)

        allocate(mesh%local_lu(FE%n_basis, FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups), &
                 mesh%local_pivots(FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups))

        !$OMP PARALLEL DO PRIVATE(mm, dir, ee, StiffnessOutflow, f, o_n, g, A, info)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:2)
            do ee = 1, mesh%n_elems
                ! SYNC FIX: Use face_normals and 0.0 threshold to match m_transport.F90 exactly
                StiffnessOutflow = -(dir(1)*mesh%elem_stiffness_x(:,:,ee) + dir(2)*mesh%elem_stiffness_y(:,:,ee))
                do f = 1, 4
                    o_n = dot_product(dir, mesh%face_normals(1:2, f, ee))
                    if (o_n > 0.0_dp) then
                        StiffnessOutflow = StiffnessOutflow + (dir(1)*mesh%face_mass_x(:,:,f,ee) + dir(2)*mesh%face_mass_y(:,:,f,ee))
                    end if
                end do

                do g = 1, n_groups
                    A = materials(mesh%material_ids(ee))%SigmaT(g) * mesh%elem_mass_matrix(:,:,ee) + StiffnessOutflow
                    
                    call dgetrf(FE%n_basis, FE%n_basis, A, FE%n_basis, mesh%local_pivots(:,ee,mm,g), info)
                    
                    if (info /= 0) then
                        write(*,*) "FATAL: LU Factorization failed (Singular Matrix) in Elem:", ee, " Angle:", mm
                        stop
                    end if

                    mesh%local_lu(:,:,ee,mm,g) = A
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Precompute_LU_Decomposition

end module m_transport_precompute
