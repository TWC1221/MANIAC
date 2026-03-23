module m_transport
    use m_constants
    use m_types
    use m_quadrature, only: t_sn_quadrature
    use m_finite_elements
    use m_material
    use m_sweep_order, only: generate_sweep_order_and_colors

    implicit none

    private
    public :: Transport_Sweep, Update_Scattering_Source_DGFEM, &
              Update_Fission_Source_DGFEM, Calculate_Total_Production_DGFEM, &
              Assemble_Fixed_Source_DGFEM, Precompute_Transport_Integrals

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(n), info
            real(8), intent(inout) :: a(lda,n), b(ldb,nrhs)
        end subroutine dgesv
    end interface

contains

subroutine Precompute_Transport_Integrals(mesh, FE, Quad, QuadBound)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad, QuadBound

        integer :: ee, q, i, j, f, i_face, j_face
        real(dp) :: elem_nodes(FE%n_basis, 2), dN_dx(FE%n_basis), dN_dy(FE%n_basis)
        real(dp) :: detJ, dV, dS
        real(dp) :: face_nodes_coords(FE%n_nodes_per_face, 2)
        real(dp) :: dx_dxi_b, dy_dxi_b

        if (.not. allocated(mesh%elem_mass_matrix)) allocate(mesh%elem_mass_matrix(FE%n_basis, FE%n_basis, mesh%n_elems))
        if (.not. allocated(mesh%elem_stiff_x))     allocate(mesh%elem_stiff_x(FE%n_basis, FE%n_basis, mesh%n_elems))
        if (.not. allocated(mesh%elem_stiff_y))     allocate(mesh%elem_stiff_y(FE%n_basis, FE%n_basis, mesh%n_elems))
        if (.not. allocated(mesh%elem_integral_v)) allocate(mesh%elem_integral_v(FE%n_basis, mesh%n_elems))
        if (.not. allocated(mesh%face_mass_matrices)) allocate(mesh%face_mass_matrices(FE%n_basis, FE%n_basis, mesh%n_faces_per_elem, mesh%n_elems))
        
        mesh%elem_mass_matrix = 0.0_dp
        mesh%elem_stiff_x     = 0.0_dp
        mesh%elem_stiff_y     = 0.0_dp
        mesh%elem_integral_v  = 0.0_dp
        mesh%face_mass_matrices = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(i, j, q, f, i_face, j_face, elem_nodes, dN_dx, dN_dy, detJ, dV, dS, &
        !$OMP&                    face_nodes_coords, dx_dxi_b, dy_dxi_b)
        do ee = 1, mesh%n_elems
            ! --- Volumetric Integration ---
            do i = 1, FE%n_basis
                elem_nodes(i, :) = mesh%nodes(mesh%elems(ee, i), :)
            end do

            do q = 1, Quad%NoPoints
                call GetMapping(FE, q, elem_nodes, dN_dx, dN_dy, detJ)
                dV = detJ * Quad%W(q)
                do i = 1, FE%n_basis
                    mesh%elem_integral_v(i, ee) = mesh%elem_integral_v(i, ee) + FE%N(q, i) * dV
                    do j = 1, FE%n_basis
                        mesh%elem_mass_matrix(i, j, ee) = mesh%elem_mass_matrix(i, j, ee) + FE%N(q, i) * FE%N(q, j) * dV
                        mesh%elem_stiff_x(i, j, ee) = mesh%elem_stiff_x(i, j, ee) + dN_dx(i) * FE%N(q, j) * dV
                        mesh%elem_stiff_y(i, j, ee) = mesh%elem_stiff_y(i, j, ee) + dN_dy(i) * FE%N(q, j) * dV
                    end do
                end do
            end do

            ! --- Face (Surface) Integration ---
            do f = 1, mesh%n_faces_per_elem
                face_nodes_coords(:,:) = mesh%nodes(mesh%elems(ee, FE%face_node_map(:,f)), :)
                
                do q = 1, QuadBound%NoPoints
                    ! Calculate surface Jacobian (length of the 1D face segment in 2D)
                    dx_dxi_b = dot_product(FE%dN_B(q, :), face_nodes_coords(:, 1))
                    dy_dxi_b = dot_product(FE%dN_B(q, :), face_nodes_coords(:, 2))
                    dS = sqrt(dx_dxi_b**2 + dy_dxi_b**2) * QuadBound%W(q)
                    
                    do i_face = 1, FE%n_nodes_per_face
                        i = FE%face_node_map(i_face, f)
                        do j_face = 1, FE%n_nodes_per_face
                            j = FE%face_node_map(j_face, f)
                            
                            mesh%face_mass_matrices(i, j, f, ee) = mesh%face_mass_matrices(i, j, f, ee) + FE%N_B(q, i_face) * FE%N_B(q, j_face) * dS
                        end do
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Precompute_Transport_Integrals

    pure subroutine Build_Volumetric_Terms(elem_id, mesh, materials, direction, group_idx, M, K)
        integer, intent(in) :: elem_id, group_idx
        type(t_mesh), intent(in) :: mesh
        type(t_material), intent(in) :: materials(:)
        real(dp), intent(in) :: direction(2)
        real(dp), intent(out) :: M(:,:), K(:,:)

        integer :: mat_id
        real(dp) :: sigma_t

        mat_id = mesh%mats(elem_id)
        sigma_t = materials(mat_id)%SigmaT(group_idx)

        M = sigma_t * mesh%elem_mass_matrix(:,:,elem_id)
        K = - (direction(1) * mesh%elem_stiff_x(:,:,elem_id) + &
               direction(2) * mesh%elem_stiff_y(:,:,elem_id))
    end subroutine Build_Volumetric_Terms

    subroutine Build_Surface_Terms(elem_id, mesh, FE, sn_quad, direction, angular_flux_vec, angular_flux_all, C, b)
        integer, intent(in) :: elem_id
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: sn_quad
        real(dp), intent(in) :: direction(2)
        real(dp), intent(in) :: angular_flux_vec(:) 
        real(dp), intent(in) :: angular_flux_all(:,:)
        real(dp), intent(out) :: C(:,:), b(:)

        integer :: f, i, neighbor_id, neighbor_face, orientation
        integer :: j_face, j_neighbor_face_node_idx, neighbor_dof_idx
        real(dp) :: omega_dot_n, upwind_val, ref_dir(2), max_dot, dprod
        integer :: mm_ref, m_iter, my_dof_idx

        C = 0.0_dp
        b = 0.0_dp

        do f = 1, mesh%n_faces_per_elem
            omega_dot_n = direction(1)*mesh%face_normals(1,f,elem_id) + &
                          direction(2)*mesh%face_normals(2,f,elem_id)
            
            if (omega_dot_n > 0.0_dp) then 
                ! OUTGOING: scale the precomputed face mass matrix
                C = C + omega_dot_n * mesh%face_mass_matrices(:,:,f,elem_id)
            else 
                ! INCOMING: Pull from neighbor
                neighbor_id = mesh%face_info(1, f, elem_id)
                if (neighbor_id > 0) then 
                    neighbor_face = mesh%face_info(2, f, elem_id)
                    orientation   = mesh%face_info(3, f, elem_id)
                    
                    ! Note: The upwind flux accumulation still requires a node loop,
                    ! but the quadrature integration is now replaced by the mass matrix.
                    do j_face = 1, FE%n_nodes_per_face
                        if (orientation == 1) then 
                            j_neighbor_face_node_idx = j_face
                        else 
                            j_neighbor_face_node_idx = FE%n_nodes_per_face - j_face + 1
                        end if
                        
                        neighbor_dof_idx = (neighbor_id - 1) * FE%n_basis + &
                                           FE%face_node_map(j_neighbor_face_node_idx, neighbor_face)
                        
                        upwind_val = -omega_dot_n * angular_flux_vec(neighbor_dof_idx)
                        
                        ! Update RHS vector b using precomputed coupling
                        do i = 1, FE%n_basis
                            ! Coupling the neighbor's node j_face to our local nodes
                            b(i) = b(i) + upwind_val * mesh%face_mass_matrices(i, FE%face_node_map(j_face,f), f, elem_id)
                        end do
                    end do
                else
                    ! Boundary Condition
                    if (mesh%face_info(4, f, elem_id) >= 103 .and. mesh%face_info(4, f, elem_id) <= 104) then
                        ! Reflective Boundary (101-104)
                        ref_dir = direction - 2.0_dp * omega_dot_n * mesh%face_normals(:, f, elem_id)
                        
                        max_dot = -1.0_dp
                        mm_ref = 1
                        do m_iter = 1, sn_quad%NoAngles
                            dprod = dot_product(ref_dir, sn_quad%Angles(m_iter, 1:2))
                            if (dprod > max_dot) then
                                max_dot = dprod
                                mm_ref = m_iter
                            end if
                        end do
                        
                        do j_face = 1, FE%n_nodes_per_face
                            my_dof_idx = (elem_id - 1) * FE%n_basis + FE%face_node_map(j_face, f)
                            upwind_val = -omega_dot_n * angular_flux_all(my_dof_idx, mm_ref)
                            
                            do i = 1, FE%n_basis
                                b(i) = b(i) + upwind_val * mesh%face_mass_matrices(i, FE%face_node_map(j_face,f), f, elem_id)
                            end do
                        end do
                    end if
                end if
            end if
        end do 
    end subroutine Build_Surface_Terms

   subroutine Transport_Sweep(mesh, FE, materials, sn_quad, &
                                     angular_flux, scalar_flux, fixed_source, & 
                                     scatter_source, fission_source, n_groups)
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_material), intent(in) :: materials(:)
        type(t_sn_quadrature), intent(in) :: sn_quad
        real(dp), intent(inout) :: angular_flux(:,:,:) 
        real(dp), intent(inout) :: scalar_flux(:,:)    
        real(dp), intent(in) :: fixed_source(:,:)    
        real(dp), intent(in) :: scatter_source(:,:)    
        real(dp), intent(in) :: fission_source(:,:)    
        integer, intent(in) :: n_groups

        integer :: mm, g, i, ie, cc, j
        real(dp) :: dir_buffer(2) 
        integer, allocatable :: sweep_order(:,:), colors(:,:), num_colors_ang(:)

        allocate(sweep_order(mesh%n_elems, sn_quad%NoAngles))
        allocate(colors(mesh%n_elems, sn_quad%NoAngles))
        allocate(num_colors_ang(sn_quad%NoAngles))
    
        do mm = 1, sn_quad%NoAngles
            dir_buffer(1) = sn_quad%Angles(mm, 1)
            dir_buffer(2) = sn_quad%Angles(mm, 2)
            call generate_sweep_order_and_colors(mesh, dir_buffer, sweep_order(:,mm), colors(:,mm), num_colors_ang(mm))
        end do

        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
        !$OMP& PRIVATE(mm, g, cc, i, ie, j)
        do mm = 1, sn_quad%NoAngles
            block
                real(dp) :: direction(2)
                real(dp) :: M(FE%n_basis, FE%n_basis), Flow_Matrix(FE%n_basis, FE%n_basis)
                real(dp) :: A(FE%n_basis, FE%n_basis)
                real(dp) :: b(FE%n_basis)
                integer  :: ipiv(FE%n_basis), info, global_indices(FE%n_basis)
                
                integer :: f, mat_id, n_incoming
                integer :: inc_face_idx(4)
                real(dp) :: inc_omega_dot_n(4)
                real(dp) :: omega_dot_n, upwind_val, ref_dir(2), max_dot, dprod
                integer :: neighbor_id, neighbor_face, orientation
                integer :: j_face, j_neighbor_face_node_idx, neighbor_dof_idx, my_dof_idx
                integer :: mm_ref, m_iter

                direction(1) = sn_quad%Angles(mm, 1)
                direction(2) = sn_quad%Angles(mm, 2)

                do cc = 1, num_colors_ang(mm)
                    do i = 1, mesh%n_elems
                        ie = sweep_order(i, mm)
                        if (colors(ie, mm) /= cc) cycle
                        
                        global_indices = (ie - 1) * FE%n_basis + [(j, j=1, FE%n_basis)]
                        mat_id = mesh%mats(ie)

                        ! Precompute Flow Matrix (K + C_outgoing) for this angle/element
                        Flow_Matrix = - (direction(1) * mesh%elem_stiff_x(:,:,ie) + &
                                         direction(2) * mesh%elem_stiff_y(:,:,ie))

                        n_incoming = 0
                        do f = 1, mesh%n_faces_per_elem
                            omega_dot_n = direction(1)*mesh%face_normals(1,f,ie) + &
                                          direction(2)*mesh%face_normals(2,f,ie)
                            
                            if (omega_dot_n > 0.0_dp) then
                                Flow_Matrix = Flow_Matrix + omega_dot_n * mesh%face_mass_matrices(:,:,f,ie)
                            else
                                n_incoming = n_incoming + 1
                                inc_face_idx(n_incoming) = f
                                inc_omega_dot_n(n_incoming) = omega_dot_n
                            end if
                        end do

                        M = mesh%elem_mass_matrix(:,:,ie)

                        do g = 1, n_groups
                            A = Flow_Matrix + materials(mat_id)%SigmaT(g) * M
                            b = fixed_source(global_indices, g) + &
                                scatter_source(global_indices, g) + &
                                fission_source(global_indices, g)

                            ! Apply Incoming fluxes
                            do j = 1, n_incoming
                                f = inc_face_idx(j)
                                omega_dot_n = inc_omega_dot_n(j)
                                neighbor_id = mesh%face_info(1, f, ie)
                                
                                if (neighbor_id > 0) then
                                    neighbor_face = mesh%face_info(2, f, ie)
                                    orientation   = mesh%face_info(3, f, ie)
                                    do j_face = 1, FE%n_nodes_per_face
                                        if (orientation == 1) then; j_neighbor_face_node_idx = j_face
                                        else; j_neighbor_face_node_idx = FE%n_nodes_per_face - j_face + 1; end if
                                        neighbor_dof_idx = (neighbor_id - 1) * FE%n_basis + FE%face_node_map(j_neighbor_face_node_idx, neighbor_face)
                                        upwind_val = -omega_dot_n * angular_flux(neighbor_dof_idx, mm, g)
                                        b = b + upwind_val * mesh%face_mass_matrices(:, FE%face_node_map(j_face,f), f, ie)
                                    end do
                                else if (mesh%face_info(4, f, ie) >= 104 .and. mesh%face_info(4, f, ie) <= 104) then
                                    ref_dir = direction - 2.0_dp * omega_dot_n * mesh%face_normals(:, f, ie)
                                    max_dot = -1.0_dp; mm_ref = 1
                                    do m_iter = 1, sn_quad%NoAngles
                                        dprod = dot_product(ref_dir, sn_quad%Angles(m_iter, 1:2))
                                        if (dprod > max_dot) then; max_dot = dprod; mm_ref = m_iter; end if
                                    end do
                                    do j_face = 1, FE%n_nodes_per_face
                                        my_dof_idx = (ie - 1) * FE%n_basis + FE%face_node_map(j_face, f)
                                        upwind_val = -omega_dot_n * angular_flux(my_dof_idx, mm_ref, g)
                                        b = b + upwind_val * mesh%face_mass_matrices(:, FE%face_node_map(j_face,f), f, ie)
                                    end do
                                end if
                            end do

                            call dgesv(FE%n_basis, 1, A, FE%n_basis, ipiv, b, FE%n_basis, info)
                            angular_flux(global_indices, mm, g) = b
                            do j = 1, FE%n_basis
                                !$OMP ATOMIC
                                scalar_flux(global_indices(j), g) = scalar_flux(global_indices(j), g) + sn_quad%W(mm) * b(j)
                            end do
                        end do
                    end do ! elements
                end do ! colors
            end block
        end do
        !$OMP END PARALLEL DO

        deallocate(sweep_order, colors, num_colors_ang)
    end subroutine Transport_Sweep

    subroutine Update_Scattering_Source_DGFEM(scat_src, scalar_flux, materials, mesh, FE, n_groups, is_adjoint)
        real(dp), intent(inout) :: scat_src(:,:)
        real(dp), intent(in)    :: scalar_flux(:,:)
        type(t_material), intent(in) :: materials(:)
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: n_groups
        logical, intent(in) :: is_adjoint

        integer :: g_to, g_from, ee, mat_id
        real(dp) :: M_phi(FE%n_basis, n_groups)
        integer :: idx_start, idx_end

        !$OMP PARALLEL DO PRIVATE(ee, mat_id, g_to, g_from, M_phi, idx_start, idx_end)
        do ee = 1, mesh%n_elems
            mat_id = mesh%mats(ee)
            idx_start = (ee - 1) * FE%n_basis + 1
            idx_end   = ee * FE%n_basis
            
            M_phi = matmul(mesh%elem_mass_matrix(:,:,ee), scalar_flux(idx_start:idx_end, :))

            do g_to = 1, n_groups
                scat_src(idx_start:idx_end, g_to) = 0.0_dp
                do g_from = 1, n_groups
                    if (.not. is_adjoint) then
                            scat_src(idx_start:idx_end, g_to) = scat_src(idx_start:idx_end, g_to) + &
                                materials(mat_id)%SigmaS(g_from, g_to) * M_phi(:, g_from)
                    else
                        ! Transpose scattering matrix for adjoint
                            scat_src(idx_start:idx_end, g_to) = scat_src(idx_start:idx_end, g_to) + &
                                materials(mat_id)%SigmaS(g_to, g_from) * M_phi(:, g_from)
                    end if
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Update_Scattering_Source_DGFEM

    subroutine Update_Fission_Source_DGFEM(fiss_src, scalar_flux, k_eff, materials, mesh, FE, n_groups, is_adjoint)
        real(dp), intent(inout)         :: fiss_src(:,:)
        real(dp), intent(in)            :: scalar_flux(:,:)
        real(dp), intent(in)            :: k_eff
        type(t_material), intent(in)    :: materials(:)
        type(t_mesh), intent(in)        :: mesh
        type(t_finite), intent(in)      :: FE
        integer, intent(in)             :: n_groups
        logical, intent(in)             :: is_adjoint

        integer                         :: g_to, g_from, ee, mat_id, i
        real(dp)                        :: local_phi(FE%n_basis), total_fiss_contrib(FE%n_basis)
        integer                         :: global_indices(FE%n_basis)

        fiss_src = 0.0_dp
        do ee = 1, mesh%n_elems
            mat_id = mesh%mats(ee)
            global_indices = (ee - 1) * FE%n_basis + [(i, i=1,FE%n_basis)]
            total_fiss_contrib = 0.0_dp
            
            do g_from = 1, n_groups
                local_phi = scalar_flux(global_indices, g_from)
                if (.not. is_adjoint) then
                    total_fiss_contrib = total_fiss_contrib + materials(mat_id)%NuSigF(g_from) * matmul(mesh%elem_mass_matrix(:,:,ee), local_phi)
                else
                    ! Adjoint sums Chi(from)
                    total_fiss_contrib = total_fiss_contrib + materials(mat_id)%Chi(g_from) * matmul(mesh%elem_mass_matrix(:,:,ee), local_phi)
                end if
            end do

            do g_to = 1, n_groups
                if (.not. is_adjoint) then
                        fiss_src(global_indices, g_to) = fiss_src(global_indices, g_to) + &
                                          (materials(mat_id)%Chi(g_to) / k_eff) * total_fiss_contrib
                else
                    ! Adjoint distributes via NuSigF(to)
                        fiss_src(global_indices, g_to) = fiss_src(global_indices, g_to) + &
                                          (materials(mat_id)%NuSigF(g_to) / k_eff) * total_fiss_contrib
                end if
            end do
        end do
    end subroutine Update_Fission_Source_DGFEM

    subroutine Calculate_Total_Production_DGFEM(total_prod, scalar_flux, materials, mesh, FE, n_groups, is_adjoint)
        type(t_material), intent(in)    :: materials(:)
        type(t_mesh), intent(in)        :: mesh
        type(t_finite), intent(in)      :: FE
        integer, intent(in)             :: n_groups
        logical, intent(in)             :: is_adjoint
        real(dp), intent(in)            :: scalar_flux(:,:)
        real(dp), intent(out)           :: total_prod

        integer                         :: g, ee, mat_id, i
        real(dp)                        :: local_phi(FE%n_basis), elem_prod
        integer                         :: global_indices(FE%n_basis)
        real(dp)                        :: adj_phi_chi_sum, nusigf_sum

        total_prod = 0.0_dp
        do ee = 1, mesh%n_elems
            mat_id = mesh%mats(ee)
            global_indices = (ee - 1) * FE%n_basis + [(i, i=1,FE%n_basis)]
            elem_prod = 0.0_dp
            
            if (.not. is_adjoint) then
                do g = 1, n_groups
                    local_phi = scalar_flux(global_indices, g)
                    elem_prod = elem_prod + materials(mat_id)%NuSigF(g) * dot_product(local_phi, mesh%elem_integral_v(:, ee))
                end do
            else
                nusigf_sum = sum(materials(mat_id)%NuSigF)
                if (nusigf_sum > 1e-12_dp) then
                    adj_phi_chi_sum = 0.0_dp
                    do g = 1, n_groups
                        local_phi = scalar_flux(global_indices, g)
                        adj_phi_chi_sum = adj_phi_chi_sum + materials(mat_id)%Chi(g) * dot_product(local_phi, mesh%elem_integral_v(:, ee))
                    end do
                    elem_prod = adj_phi_chi_sum * nusigf_sum
                end if
            end if
            total_prod = total_prod + elem_prod
        end do
    end subroutine Calculate_Total_Production_DGFEM

    subroutine Assemble_Fixed_Source_DGFEM(fixed_src, materials, mesh, FE, n_groups)
        real(dp), intent(inout) :: fixed_src(:,:)
        type(t_material), intent(in) :: materials(:)
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: n_groups

        integer :: g, ee, mat_id, i
        integer :: global_indices(FE%n_basis)

        fixed_src = 0.0_dp
        do ee = 1, mesh%n_elems
            mat_id = mesh%mats(ee)
            global_indices = (ee - 1) * FE%n_basis + [(i, i=1,FE%n_basis)]
            do g = 1, n_groups
                if (materials(mat_id)%Src(g) > 1e-12_dp) then
                    fixed_src(global_indices, g) = fixed_src(global_indices, g) + materials(mat_id)%Src(g) * mesh%elem_integral_v(:, ee)
                end if
            end do
        end do
    end subroutine Assemble_Fixed_Source_DGFEM

end module m_transport