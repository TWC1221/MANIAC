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
    end interface

contains

    subroutine InitialiseTransport(mesh, FE, Quad2D, Quad1D, QuadSn, materials, n_groups)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad2D, Quad1D
        type(t_sn_quadrature), intent(in) :: QuadSn
        type(t_material), intent(in) :: materials(:)
        integer, intent(in) :: n_groups

        ! 1. Geometric Integrals Mass and Stiffness
        call Precompute_Transport_Integrals(mesh, FE, Quad2D, Quad1D)
        
        ! 2. Angle Mapping 
        call Precompute_Reflective_Map(mesh, QuadSn)
        
        ! 3. LU Factorization
        call Precompute_All_Local_LU(mesh, FE, QuadSn, materials, n_groups)
        
    end subroutine InitialiseTransport

    subroutine Precompute_Transport_Integrals(mesh, FE, Quad2D, Quad1D)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad2D, Quad1D
        integer :: ee, q, i, j, f
        real(dp) :: nodes(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis), detJ, dV
        real(dp) :: R(FE%n_basis), dR_dxi(FE%n_basis), dR_deta(FE%n_basis)
        real(dp) :: u_face, v_face, detJ_face, dx_dt, dy_dt, nx_dS, ny_dS
        real(dp) :: u1, u2, v1, v2

        allocate(mesh%elem_mass_matrix(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 mesh%elem_stiffness_x(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%elem_stiffness_y(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%face_mass_x(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems), &
                 mesh%face_mass_y(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems), &
                 mesh%basis_integrals_vol(FE%n_basis, mesh%n_elems))
        
        mesh%elem_mass_matrix = 0.0_dp; mesh%elem_stiffness_x = 0.0_dp; mesh%elem_stiffness_y = 0.0_dp
        mesh%face_mass_x = 0.0_dp; mesh%face_mass_y = 0.0_dp
        mesh%face_normals = 0.0_dp
        mesh%basis_integrals_vol = 0.0_dp

        !$OMP PARALLEL DO &
        !$OMP PRIVATE(nodes, q, dN_dx, dN_dy, detJ, dV, i, j, f, R, u1, u2, v1, v2, &
        !$OMP         u_face, v_face, detJ_face, dR_dxi, dR_deta, dx_dt, dy_dt, nx_dS, ny_dS)
        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, :), :)
            u1 = mesh%span_range(1,ee); u2 = mesh%span_range(2,ee)
            v1 = mesh%span_range(3,ee); v2 = mesh%span_range(4,ee)
            do q = 1, Quad2D%n_points
                call GetMapping(FE, ee, mesh, q, Quad2D, u1, u2, v1, v2, nodes, dN_dx, dN_dy, detJ, R)
                dV = detJ * Quad2D%weights(q)
                mesh%elem_mass_matrix(:,:,ee)     = mesh%elem_mass_matrix(:,:,ee)     + spread(R,2,FE%n_basis) * spread(R,1,FE%n_basis) * dV
                mesh%elem_stiffness_x(:,:,ee)     = mesh%elem_stiffness_x(:,:,ee)     + spread(dN_dx,2,FE%n_basis)    * spread(R,1,FE%n_basis) * dV
                mesh%elem_stiffness_y(:,:,ee)     = mesh%elem_stiffness_y(:,:,ee)     + spread(dN_dy,2,FE%n_basis)    * spread(R,1,FE%n_basis) * dV
                mesh%basis_integrals_vol(:,ee)    = mesh%basis_integrals_vol(:,ee)    + R * dV
            end do
            do f = 1, mesh%n_faces_per_elem
                do q = 1, Quad1D%n_points
                    select case(f)
                    case(1); u_face = 0.5_dp*((u2-u1)*Quad1D%xi(q)+(u2+u1)); v_face = v1; detJ_face = 0.5_dp*(u2-u1)
                    case(2); u_face = u2; v_face = 0.5_dp*((v2-v1)*Quad1D%xi(q)+(v2+v1)); detJ_face = 0.5_dp*(v2-v1)
                    case(3); u_face = 0.5_dp*((u2-u1)*Quad1D%xi(q)+(u2+u1)); v_face = v2; detJ_face = 0.5_dp*(u2-u1)
                    case(4); u_face = u1; v_face = 0.5_dp*((v2-v1)*Quad1D%xi(q)+(v2+v1)); detJ_face = 0.5_dp*(v2-v1)
                    end select
                    
                    call EvalNURBS2D(FE, ee, mesh, u_face, v_face, R, dR_dxi, dR_deta)
                    
                    if (f == 1 .or. f == 3) then
                        dx_dt = dot_product(dR_dxi, nodes(:,1))
                        dy_dt = dot_product(dR_dxi, nodes(:,2))
                    else
                        dx_dt = dot_product(dR_deta, nodes(:,1))
                        dy_dt = dot_product(dR_deta, nodes(:,2))
                    end if

                    ! Ensure outward normal components match m_sweep_order logic
                    if (f == 1 .or. f == 2) then
                        nx_dS =  dy_dt * detJ_face * Quad1D%weights(q)
                        ny_dS = -dx_dt * detJ_face * Quad1D%weights(q)
                    else ! f == 3 .or. f == 4
                        nx_dS = -dy_dt * detJ_face * Quad1D%weights(q)
                        ny_dS =  dx_dt * detJ_face * Quad1D%weights(q)
                    end if

                    mesh%face_normals(1, f, ee) = mesh%face_normals(1, f, ee) + nx_dS
                    mesh%face_normals(2, f, ee) = mesh%face_normals(2, f, ee) + ny_dS
                    
                    do i = 1, FE%n_basis
                        do j = 1, FE%n_basis
                            mesh%face_mass_x(i,j,f,ee) = mesh%face_mass_x(i,j,f,ee) + R(i)*R(j) * nx_dS
                            mesh%face_mass_y(i,j,f,ee) = mesh%face_mass_y(i,j,f,ee) + R(i)*R(j) * ny_dS
                        end do
                    end do
                end do
                ! Normalize to unit vector
                detJ_face = sqrt(mesh%face_normals(1,f,ee)**2 + mesh%face_normals(2,f,ee)**2)
                mesh%face_normals(1, f, ee) = mesh%face_normals(1, f, ee) / detJ_face
                mesh%face_normals(2, f, ee) = mesh%face_normals(2, f, ee) / detJ_face
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
        integer :: ee, mm, g, f, info, i, ncp_ee
        real(dp) :: A(FE%n_basis, FE%n_basis), dir(2), o_n

        allocate(mesh%local_lu(FE%n_basis, FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups), &
                 mesh%local_pivots(FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups))

        !$OMP PARALLEL DO PRIVATE(mm, dir, ee, g, A, f, o_n, info, i, ncp_ee) SCHEDULE(DYNAMIC)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:2)
            do ee = 1, mesh%n_elems
                ncp_ee = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee)
                do g = 1, n_groups
                    A = materials(mesh%material_ids(ee))%SigmaT(g) * mesh%elem_mass_matrix(:,:,ee) - &
                        (dir(1)*mesh%elem_stiffness_x(:,:,ee) + dir(2)*mesh%elem_stiffness_y(:,:,ee))
                    do f = 1, mesh%n_faces_per_elem
                        o_n = dot_product(dir, mesh%face_normals(:, f, ee))
                        if (o_n > 0.0_dp) A = A + (dir(1)*mesh%face_mass_x(:,:,f,ee) + dir(2)*mesh%face_mass_y(:,:,f,ee))
                    end do

                    ! Safety: Pad unused basis functions to prevent singular matrices
                    do i = ncp_ee + 1, FE%n_basis
                        A(i, i) = 1.0_dp
                    end do

                    call dgetrf(FE%n_basis, FE%n_basis, A, FE%n_basis, mesh%local_pivots(:,ee,mm,g), info)
                    
                    if (info /= 0) then
                        write(*,*) "FATAL: LU Factorization failed (Singular Matrix) in Elem:", ee, " Angle:", mm
                        stop
                    end if

                    if (any(A /= A)) then
                        write(*,*) "FATAL: NaN detected in LU Factors. Elem:", ee, " Angle:", mm, " Group:", g
                        stop
                    end if

                    mesh%local_lu(:,:,ee,mm,g) = A
                end do
            end do
        end do
    end subroutine Precompute_All_Local_LU

end module m_transport_precompute