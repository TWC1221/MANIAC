module m_outpVTK_dfem
    use m_constants
    use m_types
    use m_basis !**
    use m_asmg 
    use m_material
    use m_quadrature
    implicit none
    private
    public :: export_dfem_vtk


contains

    subroutine export_dfem_vtk(filename, mesh, FE, QuadSn, scalar_flux, n_groups, is_SEM, is_adjoint)
        character(len=*), intent(in) :: filename
        type(t_mesh),     intent(in) :: mesh
        type(t_finite),   intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: QuadSn
        real(dp),         intent(in) :: scalar_flux(:,:) 
        integer,          intent(in) :: n_groups
        logical,          intent(in) :: is_SEM, is_adjoint

        integer :: ee, g, i, j, unit_v
        integer :: nbasis, npts_elem, ncells_elem
        integer :: n_sub_nodes, n_sub_elems
        integer :: gid, cid, basep, p
        integer :: refine_level
        integer :: ncp_active

        real(dp) :: u1, u2, v1, v2, u_val, v_val
        real(dp), allocatable :: xi_grid(:), eta_grid(:)
        real(dp), allocatable :: N_eval(:), reordered_coords(:,:), dR_dxi(:), dR_deta(:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)
        integer,  allocatable :: Cells(:,:)

        integer :: n00, n10, n11, n01
        character(len=1024) :: full_path, adj, sem

        ! --- N-Order Logic ---
        p            = FE%order
        ! For p=2 (9-node), refine_level=3. For p=3 (16-node), refine_level=4.
        refine_level = p + 1   
        nbasis       = FE%n_basis
        npts_elem    = refine_level**2
        ncells_elem  = (refine_level-1)**2
        n_sub_nodes  = mesh%n_elems * npts_elem
        n_sub_elems  = mesh%n_elems * ncells_elem

        allocate(xi_grid(refine_level), eta_grid(refine_level))
        allocate(N_eval(nbasis), reordered_coords(nbasis, 2), dR_dxi(nbasis), dR_deta(nbasis))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, n_groups))
        allocate(Cells(n_sub_elems, 4))

        ! Define interpolation coordinates in reference space [-1, 1]
        do i = 1, refine_level
            xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / real(p, dp)
            eta_grid(i) = xi_grid(i)
        end do

        ! 1. Interpolate geometry and solution to the sub-grid
        gid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*nbasis

            u1 = mesh%span_range(1, ee); u2 = mesh%span_range(2, ee)
            v1 = mesh%span_range(3, ee); v2 = mesh%span_range(4, ee)
            ncp_active = mesh%n_cp_xi(ee) * mesh%n_cp_eta(ee)

            ! Get actual physical coordinates of the high-order nodes
            do i = 1, ncp_active
                reordered_coords(i, :) = mesh%nodes(max(1, mesh%elems(ee, i)), :)
            end do
            
            ! Loop through eta (j) then xi (i) to match standard VTK ordering
            do j = 1, refine_level
                do i = 1, refine_level
                    gid = gid + 1
                    
                    ! Map reference space [-1, 1] to knot span [u1, u2] x [v1, v2]
                    u_val = 0.5_dp * ((u2 - u1) * xi_grid(i) + (u2 + u1))
                    v_val = 0.5_dp * ((v2 - v1) * eta_grid(j) + (v2 + v1))
                    call EvalNURBS2D(FE, ee, mesh, u_val, v_val, N_eval, dR_dxi, dR_deta)
                    
                    Xp(gid,1) = dot_product(N_eval(1:ncp_active), reordered_coords(1:ncp_active, 1))
                    Xp(gid,2) = dot_product(N_eval(1:ncp_active), reordered_coords(1:ncp_active, 2))
                    Xp(gid,3) = 0.0_dp
                    
                    do g = 1, n_groups
                        Up(gid, g) = dot_product(N_eval, scalar_flux(basep+1:basep+nbasis, g))
                    end do
                end do
            end do
        end do

        ! 2. Generate sub-quad connectivity (Linear quads for visualization)
        cid = 0
        do ee = 1, mesh%n_elems
            ! Start of this element's points in the global Xp array
            basep = (ee-1)*npts_elem 
            do j = 1, refine_level - 1
                do i = 1, refine_level - 1
                    ! Node indices in the gid sequence (1-based for now)
                    n00 = basep + (j-1)*refine_level + i
                    n10 = n00 + 1
                    n01 = basep + j*refine_level + i
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
        
        if (is_adjoint) adj = "_adj"
        if (is_SEM) sem = "_s"
        
        full_path = trim(filename) // trim(sem) // trim(adj) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSn%order)) // ".vtk"

        unit_v = 101
        open(unit=unit_v, file=trim(full_path), status='replace', action='write')
        write(unit_v, '(A)') "# vtk DataFile Version 3.0"
        write(unit_v, '(A)') "N-Order FE Solution"
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
        do i = 1, n_sub_elems
            write(unit_v, '(I2)') 9 ! VTK_QUAD
        end do

        ! Material & Pin IDs
        write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
        write(unit_v, '(A)') "SCALARS Material_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do i = 1, ncells_elem
                write(unit_v, '(I10)') mesh%material_ids(ee)
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

        close(unit_v)
        deallocate(xi_grid, eta_grid, N_eval, reordered_coords, Xp, Up, Cells, dR_dxi, dR_deta)

    end subroutine export_dfem_vtk

end module m_outpVTK_dfem