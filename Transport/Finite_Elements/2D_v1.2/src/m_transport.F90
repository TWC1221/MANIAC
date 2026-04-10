module m_transport
    use m_constants
    use m_types
    use m_quadrature, only: t_sn_quadrature
    use m_finite_elements
    use m_constants, only: check_nan_scalar, check_nan_array, check_nan_matrix
    use m_material
    use m_sweep_order

    implicit none

    private
    public :: Transport_Sweep, Calculate_Total_Production_DGFEM, Source_DGFEM

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(n), info
            real(8), intent(inout) :: a(lda,n), b(ldb,nrhs)
        end subroutine dgesv
    end interface

contains

    subroutine Transport_Sweep(mesh, FE, sn_quad, ang_flux, scalar_flux, total_source, sweep_order, n_groups, ref_ID)
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: sn_quad
        real(dp), intent(inout) :: ang_flux(:,:,:), scalar_flux(:,:)
        real(dp), intent(in) :: total_source(:,:)
        integer, intent(in) :: sweep_order(:,:), ref_ID(:), n_groups

        integer :: mm, ee, ie, g, f, m_ref, info
        integer :: neighbor_elem_id, i_face_node, i_basis, upwind_dof_idx, local_dof_idx
        real(dp) :: b(FE%n_basis), o_n, dir(2), neighbor_flux
        integer :: idx_start, idx_end

        ! Reset scalar flux
        scalar_flux = 0.0_dp

        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
        !$OMP& PRIVATE(mm, dir, ee, ie, g, b, f, o_n, neighbor_elem_id, i_face_node, i_basis, upwind_dof_idx, local_dof_idx, m_ref, info, neighbor_flux, idx_start, idx_end)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:2)
            
            do ee = 1, mesh%n_elems
                ie = sweep_order(ee, mm)
                call check_nan_scalar(real(ie, dp), "sweep_order(ee,mm)", "Transport_Sweep, Elem "//int_to_str(ee)//", Angle "//int_to_str(mm))
                idx_start = (ie-1)*FE%n_basis + 1
                idx_end   = ie*FE%n_basis

                do g = 1, n_groups
                    b = total_source(idx_start:idx_end, g)

                    do f = 1, mesh%n_faces_per_elem
                        o_n = dot_product(dir, mesh%face_normals(:, f, ie))
                        call check_nan_scalar(o_n, "o_n", "Transport_Sweep, Elem "//int_to_str(ie)//", Face "//int_to_str(f)//", Angle "//int_to_str(mm))
                        
                        if (o_n < 0.0_dp) then
                            neighbor_elem_id = mesh%face_connectivity(1, f, ie)
                            call check_nan_scalar(real(neighbor_elem_id, dp), "neighbor_elem_id", "Transport_Sweep, Elem "//int_to_str(ie)//", Face "//int_to_str(f))
                            if (neighbor_elem_id > 0) then 
                                call check_nan_array(mesh%face_mass_x(:, FE%face_node_map(i_face_node, f), f, ie), "face_mass_x", "Transport_Sweep, Elem "//int_to_str(ie)//", Face "//int_to_str(f))
                                call check_nan_array(mesh%face_mass_y(:, FE%face_node_map(i_face_node, f), f, ie), "face_mass_y", "Transport_Sweep, Elem "//int_to_str(ie)//", Face "//int_to_str(f))
                                do i_face_node = 1, FE%n_nodes_per_face
                                    upwind_dof_idx = mesh%upwind_idx(i_face_node, f, ie)
                                    neighbor_flux = ang_flux(upwind_dof_idx, mm, g)
                                    ! Vector update to b
                                    b = b - neighbor_flux * (dir(1)*mesh%face_mass_x(:, FE%face_node_map(i_face_node, f), f, ie) + &
                                                             dir(2)*mesh%face_mass_y(:, FE%face_node_map(i_face_node, f), f, ie))
                                end do
                            else if (any(mesh%face_connectivity(4, f, ie) == ref_ID)) then 
                                call check_nan_array(mesh%face_mass_x(:, FE%face_node_map(i_face_node, f), f, ie), "face_mass_x", "Transport_Sweep, Elem "//int_to_str(ie)//", Face "//int_to_str(f))
                                call check_nan_array(mesh%face_mass_y(:, FE%face_node_map(i_face_node, f), f, ie), "face_mass_y", "Transport_Sweep, Elem "//int_to_str(ie)//", Face "//int_to_str(f))
                                m_ref = mesh%reflect_map(mm, f, ie)
                                do i_face_node = 1, FE%n_nodes_per_face
                                    local_dof_idx = idx_start - 1 + FE%face_node_map(i_face_node, f)
                                    neighbor_flux = ang_flux(local_dof_idx, m_ref, g)
                                    b = b - neighbor_flux * (dir(1)*mesh%face_mass_x(:, FE%face_node_map(i_face_node, f), f, ie) + &
                                                             dir(2)*mesh%face_mass_y(:, FE%face_node_map(i_face_node, f), f, ie))
                                end do
                            end if
                        end if
                    end do

                    ! 3. Precomputed Solve
                    call check_nan_array(b, "b before dgetrs", "Transport_Sweep, Elem "//int_to_str(ie)//", Angle "//int_to_str(mm)//", Group "//int_to_str(g))
                    call dgetrs('N', FE%n_basis, 1, mesh%local_lu(:,:,ie,mm,g), &
                                FE%n_basis, mesh%local_pivots(:,ie,mm,g), b, FE%n_basis, info)
                    call check_nan_array(b, "b after dgetrs", "Transport_Sweep, Elem "//int_to_str(ie)//", Angle "//int_to_str(mm)//", Group "//int_to_str(g))

                    ! 4. Store
                    ang_flux(idx_start:idx_end, mm, g) = b
                    
                    ! 5. Update Scalar Flux (Vector atomic update)
                    do i_basis = 1, FE%n_basis
                        local_dof_idx = idx_start + i_basis - 1
                        !$OMP ATOMIC
                        scalar_flux(local_dof_idx, g) = scalar_flux(local_dof_idx, g) + sn_quad%weights(mm) * b(i_basis)
                        call check_nan_scalar(scalar_flux(local_dof_idx, g), "scalar_flux update", "Transport_Sweep, Elem "//int_to_str(ie)//", Angle "//int_to_str(mm)//", Group "//int_to_str(g))
                    end do
                end do
            end do 
        end do
        !$OMP END PARALLEL DO
    end subroutine Transport_Sweep

    subroutine Source_DGFEM(total_src, scalar_flux, k_eff, materials, mesh, FE, n_groups, is_adjoint, is_eigenvalue)
        real(dp), intent(inout) :: total_src(:,:)
        real(dp), intent(in)    :: scalar_flux(:,:)
        real(dp), intent(in)            :: k_eff
        type(t_material), intent(in) :: materials(:)
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        integer, intent(in) :: n_groups
        logical, intent(in) :: is_adjoint
        logical, intent(in) :: is_eigenvalue

        integer :: g_to, ee, mat_id
        real(dp) :: M_phi(FE%n_basis, n_groups), fission_prod(FE%n_basis)
        real(dp) :: scat_vec(n_groups), fiss_vec_in(n_groups), fiss_vec_out(n_groups)
        integer :: idx_start, idx_end

        total_src = 0.0_dp
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(mesh, FE, n_groups, materials, scalar_flux, total_src, k_eff, is_adjoint, is_eigenvalue) &
        !$OMP PRIVATE(ee, mat_id, g_to, M_phi, fission_prod, idx_start, idx_end, scat_vec, fiss_vec_in, fiss_vec_out)
        do ee = 1, mesh%n_elems
            mat_id    = mesh%material_ids(ee)
            idx_start = (ee - 1) * FE%n_basis + 1
            idx_end   = ee * FE%n_basis
        
            M_phi = matmul(mesh%elem_mass_matrix(:,:,ee), scalar_flux(idx_start:idx_end, :))

            if (is_eigenvalue) then
                fiss_vec_in  = merge(materials(mat_id)%NuSigF, materials(mat_id)%Chi, .not. is_adjoint)
                fission_prod = matmul(M_phi, fiss_vec_in)
                
                fiss_vec_out = merge(materials(mat_id)%Chi, materials(mat_id)%NuSigF, .not. is_adjoint)
            end if

            do g_to = 1, n_groups
                scat_vec = merge(materials(mat_id)%SigmaS(:, g_to), materials(mat_id)%SigmaS(g_to, :), .not. is_adjoint)
                total_src(idx_start:idx_end, g_to) = matmul(M_phi, scat_vec)

                if (is_eigenvalue) then
                    total_src(idx_start:idx_end, g_to) = total_src(idx_start:idx_end, g_to) + (fiss_vec_out(g_to) / k_eff) * fission_prod
                else
                    total_src(idx_start:idx_end, g_to) = total_src(idx_start:idx_end, g_to) + materials(mat_id)%Src(g_to) * mesh%basis_integrals_vol(:, ee)
                end if
                
                call check_nan_array(total_src(idx_start:idx_end, g_to), "total_src", "Source_DGFEM, Elem "//int_to_str(ee)//", Group "//int_to_str(g_to))
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Source_DGFEM

    subroutine Calculate_Total_Production_DGFEM(total_prod, scalar_flux, materials, mesh, FE, n_groups, is_adjoint)
        real(dp), intent(out)           :: total_prod
        real(dp), intent(in)            :: scalar_flux(:,:)
        type(t_material), intent(in)    :: materials(:)
        type(t_mesh), intent(in)        :: mesh
        type(t_finite), intent(in)      :: FE
        integer, intent(in)             :: n_groups
        logical, intent(in)             :: is_adjoint

        integer                         :: g, ee, mat_id
        integer                         :: idx_start, idx_end

        total_prod = 0.0_dp
        !$OMP PARALLEL DO DEFAULT(NONE)                                      &
        !$OMP SHARED(mesh, FE, n_groups, materials, scalar_flux, is_adjoint) &
        !$OMP PRIVATE(ee, mat_id, g, idx_start, idx_end)          &
        !$OMP REDUCTION(+:total_prod)
        do ee = 1, mesh%n_elems
            if (mesh%material_ids(ee) <= 0 .or. sum(abs(mesh%elem_mass_matrix(:,:,ee))) < 1.0d-15) cycle
            mat_id = mesh%material_ids(ee)
            idx_start = (ee - 1) * FE%n_basis + 1
            idx_end = ee * FE%n_basis

            do g = 1, n_groups
                call check_nan_array(scalar_flux(idx_start:idx_end, g), "scalar_flux", "Calculate_Total_Production_DGFEM, Elem "//int_to_str(ee)//", Group "//int_to_str(g))
                call check_nan_array(mesh%basis_integrals_vol(:, ee), "basis_integrals_vol", "Calculate_Total_Production_DGFEM, Elem "//int_to_str(ee))
                total_prod = total_prod + merge(materials(mat_id)%NuSigF(g),materials(mat_id)%Chi(g),.not. is_adjoint) * dot_product(scalar_flux(idx_start:idx_end, g), mesh%basis_integrals_vol(:, ee))
            end do
            call check_nan_scalar(total_prod, "total_prod (cumulative)", "Calculate_Total_Production_DGFEM, Elem "//int_to_str(ee))
        end do
        !$OMP END PARALLEL DO
    end subroutine Calculate_Total_Production_DGFEM

end module m_transport