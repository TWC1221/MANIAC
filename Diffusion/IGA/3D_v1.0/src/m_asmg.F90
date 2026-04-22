module m_asmg
    use m_constants
    use m_types
    implicit none
    private
    public :: read_asmg_mesh

contains
 
    subroutine read_asmg_mesh(filepath, mesh)
        character(len=*), intent(in) :: filepath
        type(t_mesh), intent(inout) :: mesh
        
        integer :: unit, iostatus, i, j, k
        character(len=1024) :: line
        integer :: patch_count, edge_count, p_count, poly_order
        integer :: max_cp, max_knots, pos, val1, val2, val3, ios

        open(newunit=unit, file=filepath, status='old', action='read')
        patch_count = 0; edge_count = 0; p_count = 0; poly_order = 0
        max_cp = 0; max_knots = 0
        
        do
            read(unit, '(A)', iostat=iostatus) line
            if (iostatus /= 0) exit
            line = adjustl(line)
            
            if (index(line, '$3D_Patch_Description_Start') > 0) then
                patch_count = patch_count + 1
            else if (index(line, 'PolyOrder') > 0) then
                pos = index(line, ':')
                if (pos == 0) pos = index(line, ' ')
                read(line(pos+1:), *) poly_order
            else if (index(line, '$1D_Patch_Description_Start') > 0 .or. &
                     index(line, '$2D_Patch_Description_Start') > 0) then
                edge_count = edge_count + 1
            else if (index(line, 'POINTS') == 1) then
                read(line(7:), *) p_count
            else if (index(line, 'control_points') > 0) then
                pos = index(line, ':') + 1
                ! Robustly read 1D, 2D, or 3D dimensions
                read(line(pos:), *, iostat=ios) val1, val2, val3
                if (ios == 0) then
                    k = val1 * val2 * val3
                else
                    read(line(pos:), *, iostat=ios) val1, val2
                    if (ios == 0) then
                        k = val1 * val2
                    else
                        read(line(pos:), *, iostat=ios) val1
                        k = val1
                    end if
                end if
                if (k > max_cp) max_cp = k
            else if (index(line, 'KnotVector') > 0) then
                pos = index(line, ':')
                if (pos == 0) then
                    if (index(line, 'KnotVectorXi') > 0) then
                        pos = index(line, 'KnotVectorXi') + 11
                    else if (index(line, 'KnotVectorEta') > 0) then
                        pos = index(line, 'KnotVectorEta') + 12
                    else if (index(line, 'KnotVectorZeta') > 0) then
                        pos = index(line, 'KnotVectorZeta') + 13
                    else
                        pos = index(line, 'KnotVector') + 9
                    end if
                end if
                read(line(pos+1:), *, iostat=iostatus) k
                if (iostatus == 0 .and. k > max_knots) max_knots = k
            end if
        end do
        mesh%n_elems = patch_count
        mesh%n_nodes = p_count
        mesh%n_edges = edge_count
        mesh%nloc    = max_cp
        mesh%order   = poly_order
        mesh%dim     = 3

        allocate(mesh%nodes(mesh%n_nodes, mesh%dim))
        allocate(mesh%weights(mesh%n_nodes))
        allocate(mesh%elems(mesh%n_elems, max_cp)); mesh%elems = 0
        allocate(mesh%knot_vectors_xi(mesh%n_elems, max_knots)); mesh%knot_vectors_xi = 0.0_dp
        allocate(mesh%knot_vectors_eta(mesh%n_elems, max_knots)); mesh%knot_vectors_eta = 0.0_dp
        allocate(mesh%knot_vectors_zeta(mesh%n_elems, max_knots)); mesh%knot_vectors_zeta = 0.0_dp
        allocate(mesh%edge_knots_xi(mesh%n_edges, max_knots)); mesh%edge_knots_xi = 0.0_dp
        allocate(mesh%edge_knots_eta(mesh%n_edges, max_knots)); mesh%edge_knots_eta = 0.0_dp
        allocate(mesh%edge_knots(mesh%n_edges, max_knots)); mesh%edge_knots = 0.0_dp
        allocate(mesh%mats(mesh%n_elems))
        allocate(mesh%edge_mats(mesh%n_edges))
        allocate(mesh%edges(mesh%n_edges, max_cp)); mesh%edges = 0

        allocate(mesh%n_cp_xi(mesh%n_elems), mesh%n_cp_eta(mesh%n_elems), mesh%n_cp_zeta(mesh%n_elems))
        allocate(mesh%n_knots_xi_patch(mesh%n_elems), mesh%n_knots_eta_patch(mesh%n_elems), mesh%n_knots_zeta_patch(mesh%n_elems))
        
        allocate(mesh%n_cp_xi_edge(mesh%n_edges), mesh%n_cp_eta_edge(mesh%n_edges))
        allocate(mesh%n_knots_xi_edge(mesh%n_edges), mesh%n_knots_eta_edge(mesh%n_edges))
        allocate(mesh%n_cp_edge(mesh%n_edges), mesh%n_knots_edge(mesh%n_edges))

        mesh%n_cp_xi = 0; mesh%n_cp_eta = 0; mesh%n_cp_zeta = 0; 
        mesh%n_knots_xi_patch = 0; mesh%n_knots_eta_patch = 0; mesh%n_knots_zeta_patch = 0
        mesh%n_cp_xi_edge = 0; mesh%n_cp_eta_edge = 0;
        mesh%n_knots_xi_edge = 0; mesh%n_knots_eta_edge = 0;
        mesh%n_cp_edge = 0; mesh%n_knots_edge = 0
                
        rewind(unit)
    
        j = 0 ! 2D Patch counter
        i = 0 ! 1D Patch counter
        
        do
            read(unit, '(A)', iostat=iostatus) line
            if (iostatus /= 0) exit
            line = adjustl(line)

            ! Skip comments and empty lines at the top level
            if (len_trim(line) == 0 .or. line(1:1) == '!') cycle

            if (index(line, 'POINTS') == 1) then
                do k = 1, mesh%n_nodes
                    read(unit, *) mesh%nodes(k,1), mesh%nodes(k,2), mesh%nodes(k,3), mesh%weights(k)
                end do

            else if (index(line, '$2D_Patch_Description_Start') > 0) then
                i = i + 1
                call parse_patch_block(unit, mesh%edges(i,:), mesh%edge_mats(i), n_cp=mesh%n_cp_edge(i), &
                                     out_knots_xi=mesh%edge_knots_xi(i,:), out_knots_eta=mesh%edge_knots_eta(i,:), &
                                     n_k_xi=mesh%n_knots_xi_edge(i), n_k_eta=mesh%n_knots_eta_edge(i), &
                                     n_cp_xi=mesh%n_cp_xi_edge(i), n_cp_eta=mesh%n_cp_eta_edge(i))
                ! Sync fallback for 1D
                mesh%edge_knots(i,:) = mesh%edge_knots_xi(i,:)
                mesh%n_knots_edge(i) = mesh%n_knots_xi_edge(i)

            else if (index(line, '$3D_Patch_Description_Start') > 0) then
                j = j + 1
                call parse_patch_block(unit, mesh%elems(j,:), mesh%mats(j), &
                                     out_knots_xi=mesh%knot_vectors_xi(j,:), out_knots_eta=mesh%knot_vectors_eta(j,:), out_knots_zeta=mesh%knot_vectors_zeta(j,:), &
                                     n_k_xi=mesh%n_knots_xi_patch(j), n_k_eta=mesh%n_knots_eta_patch(j), n_k_zeta=mesh%n_knots_zeta_patch(j), &
                                     n_cp_xi=mesh%n_cp_xi(j), n_cp_eta=mesh%n_cp_eta(j), n_cp_zeta=mesh%n_cp_zeta(j))
            end if
        end do
        
        close(unit)
        
        call write_mesh_to_files(mesh)

    end subroutine read_asmg_mesh

    subroutine parse_patch_block(u, out_cp, out_id2, n_cp, out_knots_xi, out_knots_eta, out_knots_zeta, n_k_xi, n_k_eta, n_k_zeta, n_cp_xi, n_cp_eta, n_cp_zeta)
        integer, intent(in) :: u
        integer, intent(inout) :: out_cp(:)
        integer, intent(out)   :: out_id2
        integer, intent(out), optional :: n_cp, n_cp_xi, n_cp_eta, n_cp_zeta
        real(dp), intent(inout), optional :: out_knots_xi(:), out_knots_eta(:), out_knots_zeta(:)
        integer, intent(out), optional :: n_k_xi, n_k_eta, n_k_zeta
        
        character(len=1024) :: l
        integer :: ios, n_val, m, pos, val1, val2, val3
        
        ! Initialize defaults
        out_id2 = -1
        
        do
            read(u, '(A)', iostat=ios) l
            if (ios /= 0) exit
            l = adjustl(l)
            
            if (index(l, '$1D_Patch_Description_End') > 0 .or. &
                index(l, '$2D_Patch_Description_End') > 0 .or. &
                index(l, '$3D_Patch_Description_End') > 0) exit
            
            if (l == '' .or. l(1:1) == '!') cycle

            if (index(l, 'Material_ID:') > 0) then
                pos = index(l, ':') + 1
                read(l(pos:), *, iostat=ios) out_id2
            else if (index(l, 'BC:') > 0) then
                pos = index(l, ':') + 1
                read(l(pos:), *, iostat=ios) out_id2
            else if (index(l, 'control_points:') > 0) then
                pos = index(l, ':') + 1
                val1 = 1; val2 = 1; val3 = 1
                read(l(pos:), *, iostat=ios) val1, val2, val3
                if (ios /= 0) then
                    read(l(pos:), *, iostat=ios) val1, val2
                    if (ios /= 0) read(l(pos:), *, iostat=ios) val1
                end if
                n_val = val1 * val2 * val3

                if (present(n_cp_xi))   n_cp_xi   = val1
                if (present(n_cp_eta))  n_cp_eta  = val2
                if (present(n_cp_zeta)) n_cp_zeta = val3
                if (present(n_cp))      n_cp = n_val
                
                ! Consume the control point IDs
                read(u, *, iostat=ios) (out_cp(m), m=1, n_val)
                if (ios == 0) out_cp(1:n_val) = out_cp(1:n_val) + 1

            else if (index(l, 'KnotVector') > 0 .and. present(out_knots_xi)) then
                pos = index(l, ':')
                if (pos == 0) pos = index(l, 'KnotVector') + 9
                read(l(pos+1:), *, iostat=ios) n_val
                if (ios == 0 .and. n_val > 0) then
                    read(u, *, iostat=ios) (out_knots_xi(m), m=1, n_val)
                    n_k_xi = n_val
                    if (present(n_k_eta)) n_k_eta = n_val
                    if (present(n_k_zeta)) n_k_zeta = n_val

                    if (present(out_knots_eta)) out_knots_eta(1:n_val) = out_knots_xi(1:n_val)
                    if (present(out_knots_zeta)) out_knots_zeta(1:n_val) = out_knots_xi(1:n_val)
                end if
            end if
        end do

    end subroutine parse_patch_block

    subroutine write_mesh_to_files(mesh)
        type(t_mesh), intent(in) :: mesh
        integer :: u, i

        ! 1. Export Nodes and Weights
        open(newunit=u, file='../output/nodes.dat', status='replace')
        write(u, '(A)') "# ID | X | Y | Z | Weight"
        do i = 1, mesh%n_nodes
            write(u, '(I8, 4F15.8)') i, mesh%nodes(i,1), mesh%nodes(i,2), mesh%nodes(i,3), mesh%weights(i)
        end do
        close(u)

        ! 2. Export Elements (Control Points)
        open(newunit=u, file='../output/elements.dat', status='replace')
        write(u, '(A)') "# Element ID | Control Point IDs..."
        do i = 1, mesh%n_elems
            write(u, '(I8, " : ", 500I8)') i, pack(mesh%elems(i,:), mesh%elems(i,:) /= 0)
        end do
        close(u)

        ! 3. Export Edges (1D Patches)
        open(newunit=u, file='../output/edges.dat', status='replace')
        write(u, '(A)') "# Edge ID | Boundary ID | Control Point IDs..."
        do i = 1, mesh%n_edges
            write(u, '(I8, " BC:", I4, " : ", 500I8)') i, mesh%edge_mats(i), &
                  pack(mesh%edges(i,:), mesh%edges(i,:) /= 0)
        end do
        close(u)

        ! 4. Export Materials and Pin IDs
        open(newunit=u, file='../output/materials.dat', status='replace')
        write(u, '(A)') "# Elem ID | Material ID"
        do i = 1, mesh%n_elems
            write(u, '(2I10)') i, mesh%mats(i)
        end do
        close(u)

        ! 5. Export Knot Vectors
        open(newunit=u, file='../output/knot_vectors.dat', status='replace')
        write(u, '(A)') "# Elem ID : KnotVectorXi | KnotVectorEta | KnotVectorZeta (3D Elements)"
        do i = 1, mesh%n_elems
            write(u, '(I8, " : ", 500F12.6, " | ", 500F12.6, " | ", 500F12.6)') i, &
                  pack(mesh%knot_vectors_xi(i,:), mesh%knot_vectors_xi(i,:) /= 0.0_dp), &
                  pack(mesh%knot_vectors_eta(i,:), mesh%knot_vectors_eta(i,:) /= 0.0_dp), &
                  pack(mesh%knot_vectors_zeta(i,:), mesh%knot_vectors_zeta(i,:) /= 0.0_dp)
        end do
        write(u, '(A)') "# Edge ID : KnotVector (1D Elements)"
        do i = 1, mesh%n_edges
            write(u, '(I8, " : ", 500F12.6)') i, &
                  pack(mesh%edge_knots(i,:), mesh%edge_knots(i,:) /= 0.0_dp)
        end do
        close(u)

        print*, mesh%n_elems, " elements, ", mesh%n_nodes, " nodes, ", mesh%n_edges, " edges read and exported to ../output/"
        print*, mesh%n_cp_xi, " control points in Xi direction (max), ", mesh%n_cp_eta, " in Eta direction (max), ", mesh%n_cp_zeta, " in Zeta direction (max)"
        print*, mesh%n_knots_xi_patch, " knots in Xi direction (max), ", mesh%n_knots_eta_patch, " in Eta direction (max), ", mesh%n_knots_zeta_patch, " in Zeta direction (max)"
        print*, mesh%n_cp_xi_edge, " control points in Xi direction for edges (max), ", mesh%n_cp_eta_edge, " in Eta direction for edges (max)"
        print*, mesh%n_knots_xi_edge, " knots in Xi direction for edges (max), ", mesh%n_knots_eta_edge, " in Eta direction for edges (max)"
        print*, mesh%n_cp_edge, " control points for edges (max), ", mesh%n_knots_edge, " knots for edges (max)"
        print*, "Polynomial order: ", mesh%order
    end subroutine write_mesh_to_files

end module m_asmg