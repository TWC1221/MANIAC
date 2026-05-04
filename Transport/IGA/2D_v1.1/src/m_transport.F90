module m_transport
    use m_constants
    use m_types
    use m_quadrature, only: t_sn_quadrature
    use m_basis
    use m_material
    use m_sweep_order

    implicit none

    private
    public :: Transport_Sweep, Calculate_Total_Production_DGFEM, Source_DGFEM

    interface
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

    subroutine Transport_Sweep(mesh, FE, sn_quad, ang_flux, scalar_flux, total_source, sweep_order, ref_ID)
        type(t_mesh), intent(in) :: mesh
        type(t_finite), intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: sn_quad
        real(dp), intent(inout) :: ang_flux(:,:,:), scalar_flux(:,:)
        real(dp), intent(in) :: total_source(:,:)
        integer, intent(in) :: sweep_order(:,:), ref_ID(:)

        integer :: mm, ee, ie, g, f, m_ref, info, k
        integer :: neighbor_elem_id, i_face_node, upwind_dof_idx, local_dof_idx, dof_idx
        real(dp) :: b(FE%n_basis), o_n, dir(2)
        real(dp) :: face_term(FE%n_basis, FE%n_basis)
        integer :: idx_start, idx_end, f_map(FE%n_nodes_per_face)
        logical :: is_inflow(mesh%n_faces_per_elem)
        real(dp) :: w_mm

        ! Reset scalar flux
        scalar_flux = 0.0_dp

        ! (The loop structure mm -> ee -> g is standard for SN sweeps)
        !$OMP PARALLEL DO PRIVATE(mm, dir, w_mm, ee, ie, idx_start, idx_end, f, o_n, is_inflow, &
        !$OMP&                    g, b, f_map, face_term, neighbor_elem_id, i_face_node, upwind_dof_idx, m_ref, local_dof_idx, info, k, dof_idx)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:2)
            
            do ee = 1, mesh%n_elems
                ie = sweep_order(ee, mm)
                
                idx_start = (ie-1)*FE%n_basis + 1
                idx_end   = ie*FE%n_basis
                
                ! 1. Identify inflow faces once for this element/angle
                do f = 1, mesh%n_faces_per_elem
                    o_n = dot_product(dir, mesh%face_normals(1:2, f, ie))
                    is_inflow(f) = (o_n < 0.0_dp)
                end do

                do g = 1, mesh%n_groups
                    b = total_source(idx_start:idx_end, g)

                    do f = 1, mesh%n_faces_per_elem
                        if (is_inflow(f)) then
                            f_map = FE%face_node_map(:, f)
                            face_term = dir(1)*mesh%face_mass_x(:,:,f,ie) + dir(2)*mesh%face_mass_y(:,:,f,ie)
                            neighbor_elem_id = mesh%face_connectivity(1, f, ie)
                            
                            if (neighbor_elem_id > 0) then 
                                do i_face_node = 1, FE%n_nodes_per_face
                                    upwind_dof_idx = mesh%upwind_idx(i_face_node, f, ie)
                                    b = b - ang_flux(upwind_dof_idx, mm, g) * face_term(:, f_map(i_face_node))
                                end do
                            else if (mesh%face_connectivity(4, f, ie) > 0 .and. any(mesh%face_connectivity(4, f, ie) == ref_ID)) then
                                m_ref = mesh%reflect_map(mm, f, ie)
                                do i_face_node = 1, FE%n_nodes_per_face
                                    local_dof_idx = idx_start - 1 + f_map(i_face_node)
                                    b = b - ang_flux(local_dof_idx, m_ref, g) * face_term(:, f_map(i_face_node))
                                end do
                            end if
                        end if
                    end do

                    ! 2. Solve using precomputed LU
                    call dgetrs('N', FE%n_basis, 1, mesh%local_lu(:,:,ie,mm,g), &
                                FE%n_basis, mesh%local_pivots(:,ie,mm,g), b, FE%n_basis, info)

                    ! Stability Fix: Prevent unphysical negative oscillations
                    where (b < 0.0_dp) b = 0.0_dp

                    ang_flux(idx_start:idx_end, mm, g) = b
                    
                    w_mm = sn_quad%weights(mm)/(4.0_dp * pi)
                    do k = 1, FE%n_basis
                        dof_idx = idx_start + k - 1
                        !$OMP ATOMIC
                        scalar_flux(dof_idx, g) = scalar_flux(dof_idx, g) + w_mm * b(k)
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
        real(dp) :: M_phi(FE%n_basis, n_groups)
        real(dp) :: fission_rate(FE%n_basis)
        integer :: idx_start, idx_end

        total_src = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(ee, mat_id, idx_start, idx_end, M_phi, fission_rate, g_to)
        do ee = 1, mesh%n_elems
            mat_id    = mesh%material_ids(ee)
            idx_start = (ee - 1) * FE%n_basis + 1
            idx_end   = ee * FE%n_basis
        
            M_phi = matmul(mesh%elem_mass_matrix(:,:,ee), scalar_flux(idx_start:idx_end, :))

            ! Scattering Source: Vectorized across all groups
            if (is_adjoint) then
                total_src(idx_start:idx_end, :) = matmul(M_phi, transpose(materials(mat_id)%SigmaS))
            else
                total_src(idx_start:idx_end, :) = matmul(M_phi, materials(mat_id)%SigmaS)
            end if

            if (is_eigenvalue) then
                fission_rate = matmul(M_phi, merge(materials(mat_id)%NuSigF, materials(mat_id)%Chi, .not. is_adjoint)) / (k_eff)
                do g_to = 1, n_groups
                    total_src(idx_start:idx_end, g_to) = total_src(idx_start:idx_end, g_to) + fission_rate * merge(materials(mat_id)%Chi(g_to), materials(mat_id)%NuSigF(g_to), .not. is_adjoint)
                end do
            else
                do g_to = 1, n_groups
                    total_src(idx_start:idx_end, g_to) = total_src(idx_start:idx_end, g_to) + materials(mat_id)%Src(g_to) * mesh%basis_integrals_vol(:, ee)
                end do
            end if

            where (total_src(idx_start:idx_end, :) < 0.0_dp) total_src(idx_start:idx_end, :) = 0.0_dp
        end do
        !$OMP END PARALLEL DO
    end subroutine Source_DGFEM

    subroutine Calculate_Total_Production_DGFEM(total_prod, scalar_flux, materials, mesh, FE, is_adjoint)
        real(dp), intent(out)           :: total_prod
        real(dp), intent(in)            :: scalar_flux(:,:)
        type(t_material), intent(in)    :: materials(:)
        type(t_mesh), intent(in)        :: mesh
        type(t_finite), intent(in)      :: FE
        logical, intent(in)             :: is_adjoint

        integer                         :: g, ee, mat_id
        integer                         :: idx_start, idx_end

        total_prod = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(ee, mat_id, idx_start, idx_end, g) REDUCTION(+:total_prod)
        do ee = 1, mesh%n_elems
            mat_id = mesh%material_ids(ee)
            idx_start = (ee - 1) * FE%n_basis + 1
            idx_end = ee * FE%n_basis

            do g = 1, mesh%n_groups
                ! Stability: Only integrate the positive part of the flux. 
                ! This prevents numerical undershoot from collapsing the denominator.
                total_prod = total_prod + merge(materials(mat_id)%NuSigF(g),materials(mat_id)%Chi(g),.not. is_adjoint) * &
                             dot_product(max(scalar_flux(idx_start:idx_end, g), 0.0_dp), mesh%basis_integrals_vol(:, ee))
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Calculate_Total_Production_DGFEM

end module m_transport