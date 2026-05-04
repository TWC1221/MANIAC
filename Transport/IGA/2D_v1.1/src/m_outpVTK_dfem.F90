module m_outpVTK_dfem
    use m_constants
    use m_types
    use m_basis
    use m_asmg, only: derive_case_nametag, int_to_str
    use m_material
    use m_quadrature
    implicit none
    private
    public :: export_dfem_vtk

contains

    subroutine export_dfem_vtk(filename, mesh, FE, Quad, QuadSN, XPETSc, is_adjoint, refine_level_in)
        character(len=*), intent(in) :: filename
        type(t_mesh),     intent(in) :: mesh
        type(t_finite),   intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad
        type(t_sn_quadrature), intent(in) :: QuadSN
        real(dp),         intent(in) :: XPETSc(:,:) 
        logical,          intent(in) :: is_adjoint
        integer,  intent(in), optional :: refine_level_in

        integer :: ee, g, ii, jj, f, unit_v
        integer :: nbasis, npts_elem, ncells_elem
        integer :: n_sub_nodes, n_sub_elems
        integer :: gid, cid, basep, p
        integer :: refine_level
        real(dp) :: u1, u2, v1, v2, xi_val, eta_val, dR_dxi(FE%n_basis), dR_deta(FE%n_basis), detJ
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), nodes_e(FE%n_basis, 2)

        real(dp), allocatable :: xi_grid(:), eta_grid(:)
        real(dp), allocatable :: N_eval(:), dR_dxi_arr(:), dR_deta_arr(:), mesh_coords(:,:)
        real(dp), allocatable :: Xp(:,:), Up(:,:), Jp(:)
        integer,  allocatable :: Cells(:,:)

        integer :: n00, n10, n11, n01
        character(len=1024) :: full_path, mesh_path, adj, sem
        integer :: slash_idx
        integer :: pts_per_side, gid_mesh, n_mesh_pts

        call execute_command_line("mkdir -p " // trim(filename))

        ! --- N-Order Logic ---
        p            = FE%order
        if (present(refine_level_in)) then
            refine_level = refine_level_in
        else
            refine_level = p + 1
        end if
        nbasis       = FE%n_basis
        npts_elem    = refine_level**2
        ncells_elem  = (refine_level-1)**2
        n_sub_nodes  = mesh%n_elems * npts_elem
        n_sub_elems  = mesh%n_elems * ncells_elem

        allocate(xi_grid(refine_level), eta_grid(refine_level))
        allocate(N_eval(nbasis), dR_dxi_arr(nbasis), dR_deta_arr(nbasis))
        allocate(mesh_coords(nbasis, size(mesh%nodes, 2)))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, mesh%n_groups))
        allocate(Cells(n_sub_elems, 4), Jp(n_sub_elems))

        ! Define interpolation coordinates in reference space [-1, 1]
        do ii = 1, refine_level
            xi_grid(ii)  = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level - 1, dp)
            eta_grid(ii) = xi_grid(ii)
        end do

        ! 1. Interpolate geometry and solution to the sub-grid
        gid = 0
        cid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*nbasis
            u1 = mesh%elem_u_min(ee); u2 = mesh%elem_u_max(ee)
            v1 = mesh%elem_v_min(ee); v2 = mesh%elem_v_max(ee)
            ! Get actual physical coordinates of the high-order nodes
            do ii = 1, nbasis
                mesh_coords(ii, :) = mesh%nodes(mesh%elems(ee, ii), :)
            end do
            
            do jj = 1, refine_level
                do ii = 1, refine_level
                    gid = gid + 1
                    xi_val  = 0.5_dp * ((u2 - u1) * xi_grid(ii) + (u2 + u1))
                    eta_val = 0.5_dp * ((v2 - v1) * eta_grid(jj) + (v2 + v1))
                    call EvalNURBS2D(FE, ee, mesh, xi_val, eta_val, N_eval, dR_dxi_arr, dR_deta_arr)
                    
                    Xp(gid,1) = dot_product(N_eval, mesh_coords(:, 1))
                    Xp(gid,2) = dot_product(N_eval, mesh_coords(:, 2))
                    if (size(mesh_coords, 2) == 3) then
                        Xp(gid,3) = dot_product(N_eval, mesh_coords(:, 3))
                    else
                        Xp(gid,3) = 0.0_dp
                    end if
                    
                    do g = 1, mesh%n_groups
                        Up(gid, g) = dot_product(N_eval, XPETSc(basep+1:basep+nbasis, g))
                    end do
                end do
            end do
            
            ! Store Jacobian at the center of the element for debugging
            nodes_e = mesh%nodes(mesh%elems(ee, 1:nbasis), :)
            call GetMapping2D(FE, ee, mesh, 1, Quad, u1, u2, v1, v2, nodes_e, dN_dx, dN_dy, detJ, N_eval, xi_custom=0.0_dp, eta_custom=0.0_dp)
            do ii = 1, ncells_elem
                Jp(cid + ii) = detJ
            end do
            cid = cid + ncells_elem
        end do

        ! 2. Generate sub-quad connectivity (Linear quads for visualization)
        cid = 0
        do ee = 1, mesh%n_elems
            ! Start of this element's points in the global Xp array
            basep = (ee-1)*npts_elem 
            do jj = 1, refine_level - 1
                do ii = 1, refine_level - 1
                    ! Node indices in the gid sequence (1-based for now)
                    n00 = basep + (jj-1)*refine_level + ii
                    n10 = n00 + 1
                    n01 = basep + jj*refine_level + ii
                    n11 = n01 + 1
                    
                    cid = cid + 1
                    ! Store as counter-clockwise: bottom-left, bottom-right, top-right, top-left
                    Cells(cid, :) = [n00, n10, n11, n01]
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

        write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*5
        do cid = 1, n_sub_elems
            write(unit_v, '(5I10)') 4, Cells(cid,1)-1, Cells(cid,2)-1, &
                                       Cells(cid,3)-1, Cells(cid,4)-1
        end do

        write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
        do ii = 1, n_sub_elems
            write(unit_v, '(I2)') 9 ! VTK_QUAD
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
            ! Find the first non-internal/non-zero boundary ID on any face
            f = 1
            do gid = 1, 4
                if (mesh%face_connectivity(1, gid, ee) == -1) f = max(f, mesh%face_connectivity(4, gid, ee))
            end do
            do ii = 1, ncells_elem
                write(unit_v, '(I10)') f
            end do
        end do

        ! write(unit_v, '(A)') "SCALARS Jacobian_Det double 1"
        ! write(unit_v, '(A)') "LOOKUP_TABLE default"
        ! do ii = 1, n_sub_elems
        !     write(unit_v, '(F18.10)') Jp(ii)
        ! end do

        ! Solution Data
        write(unit_v, '(A, I10)') "POINT_DATA ", n_sub_nodes
        do g = 1, mesh%n_groups
            write(unit_v, '(A, I0)') "SCALARS Flux_Group_", g
            write(unit_v, '(A)') "double 1"
            write(unit_v, '(A)') "LOOKUP_TABLE default"
            do gid = 1, n_sub_nodes
                write(unit_v, '(F18.10)') Up(gid,g)
            end do
        end do

        close(unit_v)

        ! --- Write Knot Span Outlines (The "Portrayal" of Elements) ---
        ! This creates a separate wireframe mesh so boundaries are clearly visible
        pts_per_side = max(refine_level, 10) 
        n_mesh_pts   = mesh%n_elems * pts_per_side * 4
        
        open(unit=unit_v, file=trim(mesh_path), status='replace', action='write')
        write(unit_v, '(A)') "# vtk DataFile Version 3.0"
        write(unit_v, '(A)') "IGA Knot Span Boundaries"
        write(unit_v, '(A)') "ASCII"
        write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"
        
        write(unit_v, '(A, I0, A)') "POINTS ", n_mesh_pts, " double"
        do ee = 1, mesh%n_elems
            u1 = mesh%elem_u_min(ee); u2 = mesh%elem_u_max(ee)
            v1 = mesh%elem_v_min(ee); v2 = mesh%elem_v_max(ee)
            mesh_coords(1:nbasis, :) = mesh%nodes(mesh%elems(ee, 1:nbasis), :)
            
            ! Loop over 4 sides of the knot span
            do f = 1, 4
                do ii = 1, pts_per_side
                    xi_val = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(pts_per_side-1, dp)
                    select case(f)
                        case(1); xi_val =  xi_val; eta_val = -1.0_dp
                        case(2); eta_val = xi_val; xi_val =   1.0_dp
                        case(3); xi_val = -xi_val; eta_val =  1.0_dp
                        case(4); eta_val = -xi_val; xi_val = -1.0_dp
                    end select
                    
                    ! Map to parametric space
                    xi_val  = 0.5_dp * ((u2 - u1) * xi_val + (u2 + u1))
                    eta_val = 0.5_dp * ((v2 - v1) * eta_val + (v2 + v1))
                    
                    call EvalNURBS2D(FE, ee, mesh, xi_val, eta_val, N_eval, dR_dxi_arr, dR_deta_arr)
                    write(unit_v, '(3F18.10)') dot_product(N_eval, mesh_coords(:, 1)), &
                                               dot_product(N_eval, mesh_coords(:, 2)), 0.0_dp
                end do
            end do
        end do
        
        write(unit_v, '(A, 2I10)') "CELLS ", mesh%n_elems, mesh%n_elems * (pts_per_side * 4 + 2)
        gid_mesh = 0
        do ee = 1, mesh%n_elems
            write(unit_v, '(I10)', advance='no') pts_per_side * 4 + 1
            do ii = 1, pts_per_side * 4
                write(unit_v, '(I10)', advance='no') gid_mesh
                gid_mesh = gid_mesh + 1
            end do
            ! Close the loop by connecting back to the start of this element's points
            write(unit_v, '(I10)') gid_mesh - (pts_per_side * 4)
        end do
        
        write(unit_v, '(A, I10)') "CELL_TYPES ", mesh%n_elems
        do ee = 1, mesh%n_elems
            write(unit_v, '(I2)') 4 ! VTK_POLY_LINE
        end do
        
        write(unit_v, '(A, I10)') "CELL_DATA ", mesh%n_elems
        write(unit_v, '(A)') "SCALARS Element_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            write(unit_v, '(I10)') ee
        end do
        
        write(unit_v, '(A)') "SCALARS Patch_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            write(unit_v, '(I10)') mesh%elem_patch_id(ee)
        end do

        write(unit_v, '(A)') "SCALARS Material_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            write(unit_v, '(I10)') mesh%material_ids(ee)
        end do
        
        close(unit_v)

        deallocate(xi_grid, eta_grid, N_eval, dR_dxi_arr, dR_deta_arr, mesh_coords, Xp, Up, Cells)
    end subroutine export_dfem_vtk

end module m_outpVTK_dfem