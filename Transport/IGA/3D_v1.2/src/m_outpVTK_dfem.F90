module m_outpVTK_dfem
    use m_constants
    use m_types
    use m_basis
    use m_asmg, only: derive_case_nametag, int_to_str
    use m_material
    use m_quadrature
    use iso_fortran_env, only: error_unit
    implicit none
    private
    public :: export_dfem_vtk

contains

    subroutine export_dfem_vtk(filename, mesh, FE, Quad, QuadSN, XPETSc, is_adjoint, refine_level_in, ang_flux_in)
        character(len=*), intent(in) :: filename
        type(t_mesh),     intent(in) :: mesh
        type(t_finite),   intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad
        type(t_sn_quadrature), intent(in) :: QuadSN
        real(dp),         intent(in) :: XPETSc(:,:) 
        logical,          intent(in) :: is_adjoint
        integer,  intent(in), optional :: refine_level_in
        real(dp), intent(in), optional :: ang_flux_in(:,:,:)

        integer :: ee, g, ii, jj, kk, f, unit_v
        integer :: nbasis, npts_elem, ncells_elem, mm
        integer :: n_sub_nodes, n_sub_elems, n_angles
        integer :: gid, cid, basep, p, n_groups, current_angle_idx
        real(dp) :: rand_val
        integer, parameter :: export_angle_limit = 5 ! Limit to export only the first 5 angles
        integer :: refine_level
        real(dp) :: u1, u2, v1, v2, w1, w2, xi_val, eta_val, zeta_val

        real(dp), allocatable :: xi_grid(:), eta_grid(:), zeta_grid(:)
        real(dp), allocatable :: N_eval(:), dR_dxi_arr(:), dR_deta_arr(:), dR_dzeta_arr(:), mesh_coords(:,:)
        real(dp), allocatable :: Xp(:,:), Up(:,:), Jp(:), AngUp(:,:,:)
        integer,  allocatable :: Cells(:,:)

        integer, allocatable :: angles_to_export_indices(:)
        integer :: n000, n100, n110, n010, n001, n101, n111, n011
        character(len=1024) :: full_path, mesh_path, adj, sem
        integer :: slash_idx, j_rand, n_angles_to_export

        call execute_command_line("mkdir -p " // trim(filename))

        ! --- N-Order Logic ---
        p            = FE%order
        if (present(refine_level_in)) then
            refine_level = refine_level_in
        else
            refine_level = p + 1
        end if
        nbasis       = FE%n_basis
        npts_elem    = refine_level**3
        ncells_elem  = (refine_level-1)**3
        n_sub_nodes  = mesh%n_elems * npts_elem
        n_sub_elems  = mesh%n_elems * ncells_elem
        n_groups     = mesh%n_groups
        n_angles     = QuadSN%n_angles

        allocate(xi_grid(refine_level), eta_grid(refine_level), zeta_grid(refine_level))
        allocate(N_eval(nbasis), dR_dxi_arr(nbasis), dR_deta_arr(nbasis), dR_dzeta_arr(nbasis))
        allocate(mesh_coords(nbasis, size(mesh%nodes, 2)))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, n_groups))
        allocate(Cells(n_sub_elems, 8), Jp(n_sub_elems))

        if (present(ang_flux_in)) then
            allocate(AngUp(n_sub_nodes, n_angles, n_groups))

            ! --- Select random angles for export ---
            ! Seed the random number generator (for different results each run)
            call random_seed()

            ! Determine how many angles to export (min of limit and actual angles)
            n_angles_to_export = min(export_angle_limit, n_angles)
            allocate(angles_to_export_indices(n_angles_to_export))

            if (n_angles_to_export == n_angles) then
                ! If exporting all available angles, just use them in order
                do ii = 1, n_angles_to_export
                    angles_to_export_indices(ii) = ii
                end do
            else
                ! Select unique random angles
                j_rand = 0
                do while (j_rand < n_angles_to_export)
                    call random_number(rand_val)
                    ! Generate a random index between 1 and n_angles
                    current_angle_idx = floor(rand_val * n_angles) + 1
                    
                    ! Check for uniqueness
                    if (.not. any(angles_to_export_indices(1:j_rand) == current_angle_idx)) then
                        j_rand = j_rand + 1
                        angles_to_export_indices(j_rand) = current_angle_idx
                    end if
                end do
            end if

        end if

        ! Define interpolation coordinates in reference space [-1, 1]
        do ii = 1, refine_level
            xi_grid(ii)  = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level - 1, dp)
            eta_grid(ii) = xi_grid(ii)
            zeta_grid(ii) = xi_grid(ii)
        end do

        ! 1. Interpolate geometry and solution to the sub-grid
        gid = 0
        cid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*nbasis
            u1 = mesh%elem_u_min(ee); u2 = mesh%elem_u_max(ee)
            v1 = mesh%elem_v_min(ee); v2 = mesh%elem_v_max(ee)
            w1 = mesh%elem_w_min(ee); w2 = mesh%elem_w_max(ee)

            ! Get actual physical coordinates of the high-order nodes
            do ii = 1, nbasis
                mesh_coords(ii, :) = mesh%nodes(mesh%elems(ee, ii), :)
            end do
            
            do kk = 1, refine_level
                do jj = 1, refine_level
                    do ii = 1, refine_level
                        gid = gid + 1
                        xi_val   = 0.5_dp * ((u2 - u1) * xi_grid(ii) + (u2 + u1))
                        eta_val  = 0.5_dp * ((v2 - v1) * eta_grid(jj) + (v2 + v1))
                        zeta_val = 0.5_dp * ((w2 - w1) * zeta_grid(kk) + (w2 + w1))
                        
                        call EvalNURBS3D(FE, ee, mesh, xi_val, eta_val, zeta_val, &
                                         N_eval, dR_dxi_arr, dR_deta_arr, dR_dzeta_arr)
                        
                        Xp(gid,1) = dot_product(N_eval, mesh_coords(:, 1))
                        Xp(gid,2) = dot_product(N_eval, mesh_coords(:, 2))
                        Xp(gid,3) = dot_product(N_eval, mesh_coords(:, 3))
                        
                        do g = 1, n_groups
                            Up(gid, g) = dot_product(N_eval, XPETSc(basep+1:basep+nbasis, g))
                        end do

                        if (present(ang_flux_in)) then
                            ! Iterate over the selected random angles
                            do mm = 1, n_angles_to_export
                                do g = 1, n_groups
                                    AngUp(gid, angles_to_export_indices(mm), g) = dot_product(N_eval, ang_flux_in(basep+1:basep+nbasis, angles_to_export_indices(mm), g))
                                end do
                            end do
                        end if
                    end do
                end do
            end do
        end do

        ! 2. Generate sub-hexahedron connectivity
        cid = 0
        do ee = 1, mesh%n_elems
            ! Start of this element's points in the global Xp array
            basep = (ee-1)*npts_elem 
            do kk = 1, refine_level - 1
                do jj = 1, refine_level - 1
                    do ii = 1, refine_level - 1
                        n000 = basep + (kk-1)*refine_level**2 + (jj-1)*refine_level + ii
                        n100 = n000 + 1
                        n110 = n000 + refine_level + 1
                        n010 = n000 + refine_level
                        n001 = n000 + refine_level**2
                        n101 = n001 + 1
                        n111 = n001 + refine_level + 1
                        n011 = n001 + refine_level
                        
                        cid = cid + 1
                        Cells(cid, :) = [n000, n100, n110, n010, n001, n101, n111, n011]
                    end do
                end do
            end do
        end do

        ! 4. File Writing
        adj = ""
        sem = "_f"
        slash_idx = index(filename, '/', back=.true.)
        if (slash_idx == 0) slash_idx = index(filename, '\', back=.true.)
        
        if (is_adjoint) adj = "_adj"
        
        ! Field File Path
        if (slash_idx > 0) then
            full_path = trim(filename) // "/" // trim(filename(slash_idx+1:)) // trim(sem) // trim(adj) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSN%order)) // ".vtk"
            mesh_path = trim(filename) // "/" // trim(filename(slash_idx+1:)) // "_mesh" // ".vtk"
        else
            full_path = trim(filename) // "/" // trim(filename) // trim(sem) // trim(adj) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSN%order)) // ".vtk"
            mesh_path = trim(filename) // "/" // trim(filename) // "_mesh" // ".vtk"
        end if

        ! --- Write Field Data ---
        unit_v = 101
        open(unit=unit_v, file=trim(full_path), status='replace', action='write')
        write(unit_v, '(A)') "# vtk DataFile Version 3.0"
        write(unit_v, '(A)') "IGA DGFEM Solution Field"
        write(unit_v, '(A)') "ASCII"
        write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

        write(unit_v, '(A, I0, A)') "POINTS ", n_sub_nodes, " double"
        do gid = 1, n_sub_nodes
            write(unit_v, '(3F18.10)') Xp(gid,:)
        end do

        write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*9
        do cid = 1, n_sub_elems
            write(unit_v, '(9I10)') 8, Cells(cid,1)-1, Cells(cid,2)-1, Cells(cid,3)-1, Cells(cid,4)-1, &
                                       Cells(cid,5)-1, Cells(cid,6)-1, Cells(cid,7)-1, Cells(cid,8)-1
        end do

        write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
        do ii = 1, n_sub_elems
            write(unit_v, '(I2)') 12 ! VTK_HEXAHEDRON
        end do

        ! Material & Pin IDs
        write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
        
        write(unit_v, '(A)') "SCALARS Material_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, ncells_elem
                write(unit_v, '(I10)') mesh%material_ids(ee)
            end do
        end do

        write(unit_v, '(A)') "SCALARS Patch_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, ncells_elem
                write(unit_v, '(I10)') mesh%elem_patch_id(ee)
            end do
        end do

        ! Consolidate BC Faces into a single Boundary_ID field
        write(unit_v, '(A)') "SCALARS Boundary_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            ! Identify if this element sits on a boundary
            f = 0 ! 0 = Internal
            do gid = 1, mesh%n_faces_per_elem
                if (mesh%face_connectivity(1, gid, ee) == -1) then
                    f = max(f, mesh%face_connectivity(4, gid, ee))
                end if
            end do
            
            do ii = 1, ncells_elem
                write(unit_v, '(I10)') f
            end do
        end do

        ! Solution Data
        write(unit_v, '(A, I10)') "POINT_DATA ", n_sub_nodes
        do g = 1, n_groups
            write(unit_v, '(A, I0)') "SCALARS Flux_Group_", g
            write(unit_v, '(A)') "double 1"
            write(unit_v, '(A)') "LOOKUP_TABLE default"
            do gid = 1, n_sub_nodes
                write(unit_v, '(F18.10)') Up(gid,g)
            end do
        end do

        if (present(ang_flux_in)) then
            ! Loop over the selected random angles
            do mm = 1, n_angles_to_export
                do g = 1, n_groups
                    write(unit_v, '(A, I0, A, I0)') "SCALARS Angular_Flux_G", g, "_Ang", angles_to_export_indices(mm)
                    write(unit_v, '(A)') "double 1"
                    write(unit_v, '(A)') "LOOKUP_TABLE default"
                    do gid = 1, n_sub_nodes
                        write(unit_v, '(F18.10)') AngUp(gid, angles_to_export_indices(mm), g)
                    end do
                end do
            end do
        end if

        close(unit_v)

        deallocate(xi_grid, eta_grid, zeta_grid, N_eval, dR_dxi_arr, dR_deta_arr, dR_dzeta_arr, mesh_coords, Xp, Up, Cells, angles_to_export_indices)
        if (allocated(AngUp)) deallocate(AngUp)
    end subroutine export_dfem_vtk

end module m_outpVTK_dfem