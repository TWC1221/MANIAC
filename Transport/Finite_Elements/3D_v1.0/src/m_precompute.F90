module m_transport_precompute
    use m_types
    use m_quadrature
    use m_finite_elements
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

    subroutine InitialiseTransport(mesh, FE, Quad3D, Quad2D, sn_quad, materials, n_groups)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad3D, Quad2D
        type(t_sn_quadrature), intent(in) :: sn_quad
        type(t_material), intent(in) :: materials(:)
        integer, intent(in) :: n_groups

        call timer_start('Precomputation')
        ! 1. Geometric Integrals Mass and Stiffness
        call timer_start('EVAL - Transport Integrals')
        call Precompute_Transport_Integrals(mesh, FE, Quad3D, Quad2D)
        call timer_stop('EVAL - Transport Integrals')
        ! 2. Angle Mapping 
        call timer_start('EVAL - Reflective BC Mapping')
        call Precompute_Reflective_Map(mesh, sn_quad)
        call timer_stop('EVAL - Reflective BC Mapping')
        ! 3. LU Factorization
        call timer_start('EVAL - LU Decomposition')
        call Precompute_All_Local_LU(mesh, FE, sn_quad, materials, n_groups)
        call timer_stop('EVAL - LU Decomposition')

        call timer_stop('Precomputation')
        
    end subroutine InitialiseTransport

    subroutine Precompute_Transport_Integrals(mesh, FE, Quad3D, Quad2D)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad3D, Quad2D
        integer :: ee, q, i, j, f, i_f, j_f
        real(dp) :: nodes(FE%n_basis, 3), dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis), detJ, dV
        real(dp) :: coords(FE%n_nodes_per_face, 3)
        real(dp) :: t1(3), t2(3), n_vec(3), elem_centroid(3), q_phys(3), vec_diff(3)
        real(dp) :: w

        allocate(mesh%elem_mass_matrix(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 mesh%elem_stiffness_x(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%elem_stiffness_y(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%elem_stiffness_z(FE%n_basis, FE%n_basis, mesh%n_elems),     &
                 mesh%face_mass_x(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems), &
                 mesh%face_mass_y(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems), &
                 mesh%face_mass_z(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems), &
                 mesh%basis_integrals_vol(FE%n_basis, mesh%n_elems))
        
        mesh%elem_mass_matrix = 0.0_dp
        mesh%elem_stiffness_x = 0.0_dp; mesh%elem_stiffness_y = 0.0_dp; mesh%elem_stiffness_z = 0.0_dp
        mesh%face_mass_x = 0.0_dp; mesh%face_mass_y = 0.0_dp; mesh%face_mass_z = 0.0_dp
        mesh%basis_integrals_vol = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(nodes, q, dN_dx, dN_dy, dN_dz, detJ, dV, i, j, f, coords, t1, t2, n_vec, i_f, j_f, elem_centroid, q_phys, vec_diff, w)
        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, :), :)
            
            ! Calculate element centroid for outward normal check
            elem_centroid = sum(nodes, dim=1) / real(FE%n_basis, dp)
            
            do q = 1, Quad3D%n_points
                call GetMapping(FE, q, nodes, dN_dx, dN_dy, dN_dz, detJ)
                dV = detJ * Quad3D%weights(q)
                mesh%elem_mass_matrix(:,:,ee)     = mesh%elem_mass_matrix(:,:,ee)     + spread(FE%basis_at_quad(q,:),2,FE%n_basis) * spread(FE%basis_at_quad(q,:),1,FE%n_basis) * dV
                mesh%elem_stiffness_x(:,:,ee)     = mesh%elem_stiffness_x(:,:,ee)     + spread(dN_dx,2,FE%n_basis)    * spread(FE%basis_at_quad(q,:),1,FE%n_basis) * dV
                mesh%elem_stiffness_y(:,:,ee)     = mesh%elem_stiffness_y(:,:,ee)     + spread(dN_dy,2,FE%n_basis)    * spread(FE%basis_at_quad(q,:),1,FE%n_basis) * dV
                mesh%elem_stiffness_z(:,:,ee)     = mesh%elem_stiffness_z(:,:,ee)     + spread(dN_dz,2,FE%n_basis)    * spread(FE%basis_at_quad(q,:),1,FE%n_basis) * dV
                mesh%basis_integrals_vol(:,ee)    = mesh%basis_integrals_vol(:,ee)    + FE%basis_at_quad(q,:) * dV
            end do
            do f = 1, mesh%n_faces_per_elem
                coords = mesh%nodes(mesh%elems(ee, FE%face_node_map(:,f)), :)
                do q = 1, Quad2D%n_points
                    ! Calculate tangents and physical position on face
                    t1(1) = dot_product(coords(:,1), FE%dbasis_bound_dxi(q,:)); t2(1) = dot_product(coords(:,1), FE%dbasis_bound_deta(q,:))
                    t1(2) = dot_product(coords(:,2), FE%dbasis_bound_dxi(q,:)); t2(2) = dot_product(coords(:,2), FE%dbasis_bound_deta(q,:))
                    t1(3) = dot_product(coords(:,3), FE%dbasis_bound_dxi(q,:)); t2(3) = dot_product(coords(:,3), FE%dbasis_bound_deta(q,:))
                    
                    q_phys(1) = dot_product(coords(:,1), FE%basis_at_bound_quad(q,:))
                    q_phys(2) = dot_product(coords(:,2), FE%basis_at_bound_quad(q,:))
                    q_phys(3) = dot_product(coords(:,3), FE%basis_at_bound_quad(q,:))

                    ! n_vec here is unnormalized normal * dS (Jacobian included)
                    n_vec(1) = t1(2)*t2(3) - t1(3)*t2(2)
                    n_vec(2) = t1(3)*t2(1) - t1(1)*t2(3)
                    n_vec(3) = t1(1)*t2(2) - t1(2)*t2(1)

                    ! Ensure outward pointing
                    vec_diff = q_phys - elem_centroid
                    if (dot_product(n_vec, vec_diff) < 0.0_dp) n_vec = -n_vec
                    
                    w = Quad2D%weights(q)

                    do i_f = 1, FE%n_nodes_per_face
                        do j_f = 1, FE%n_nodes_per_face
                            i = FE%face_node_map(i_f, f); j = FE%face_node_map(j_f, f)
                            mesh%face_mass_x(i,j,f,ee) = mesh%face_mass_x(i,j,f,ee) + FE%basis_at_bound_quad(q,i_f)*FE%basis_at_bound_quad(q,j_f) * n_vec(1) * w
                            mesh%face_mass_y(i,j,f,ee) = mesh%face_mass_y(i,j,f,ee) + FE%basis_at_bound_quad(q,i_f)*FE%basis_at_bound_quad(q,j_f) * n_vec(2) * w
                            mesh%face_mass_z(i,j,f,ee) = mesh%face_mass_z(i,j,f,ee) + FE%basis_at_bound_quad(q,i_f)*FE%basis_at_bound_quad(q,j_f) * n_vec(3) * w
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
                normal = mesh%face_normals(:, f, ee)
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
        real(dp) :: A(FE%n_basis, FE%n_basis), dir(3), o_n

        allocate(mesh%local_lu(FE%n_basis, FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups), &
                 mesh%local_pivots(FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups))

        !$OMP PARALLEL DO PRIVATE(mm, dir, ee, g, A, f, o_n, info) SCHEDULE(DYNAMIC)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, :)
            do ee = 1, mesh%n_elems
                do g = 1, n_groups
                    A = materials(mesh%material_ids(ee))%SigmaT(g) * mesh%elem_mass_matrix(:,:,ee) - &
                        (dir(1)*mesh%elem_stiffness_x(:,:,ee) + dir(2)*mesh%elem_stiffness_y(:,:,ee) + dir(3)*mesh%elem_stiffness_z(:,:,ee))
                    do f = 1, mesh%n_faces_per_elem
                        o_n = dot_product(dir, mesh%face_normals(:, f, ee))
                        if (o_n > 0.0_dp) A = A + (dir(1)*mesh%face_mass_x(:,:,f,ee) + dir(2)*mesh%face_mass_y(:,:,f,ee) + dir(3)*mesh%face_mass_z(:,:,f,ee))
                    end do
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