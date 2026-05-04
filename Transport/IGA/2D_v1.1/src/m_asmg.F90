module m_asmg
    use m_constants
    use m_types
    implicit none
    private
    public :: read_asmg_mesh, derive_case_nametag, int_to_str

contains

    pure function int_to_str(ii) result(res)
        integer, intent(in) :: ii
        character(len=12) :: res
        write(res, '(I0)') ii
        res = adjustl(res)
    end function int_to_str

    pure function derive_case_nametag(filename) result(tag)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: tag
        integer :: pos, dot_pos

        pos = index(filename, '/', back=.true.)
        if (pos == 0) pos = index(filename, '\', back=.true.)

        tag = filename(pos+1:)
        dot_pos = index(tag, '.', back=.true.)
        
        if (dot_pos > 0) then
            tag = tag(1:dot_pos-1) // ".vtk"
        else
            tag = tag // ".vtk"
        end if
    end function derive_case_nametag
 
    subroutine read_asmg_mesh(filepath, mesh)
        character(len=*), intent(in) :: filepath
        type(t_mesh), intent(inout) :: mesh

        integer :: unit, iostatus, ii, jj, k, p_idx, span_u, span_v, cp_idx
        real(dp) :: dummy_coord
        character(len=1024) :: line
        integer :: patch_count, edge_count, n_total_elems
        integer :: max_cp, max_knots, pos, p, q

        open(newunit=unit, file=filepath, status='old', action='read')
        mesh%order     = 1 ! Default polynomial order
        mesh%n_nodes  = 0
        patch_count = 0; edge_count = 0;
        max_cp = 0; max_knots = 0
        
        do
            read(unit, '(A)', iostat=iostatus) line
            if (iostatus /= 0) exit
            line = adjustl(line)
            
            if (index(line, '$2D_Patch_Description_Start') > 0) then
                patch_count = patch_count + 1
            else if (index(line, 'PolyOrder') > 0) then
                pos = index(line, 'PolyOrder') + 9
                read(line(pos+1:), *) mesh%order
            else if (index(line, 'Groups') > 0) then
                pos = index(line, 'Groups') + 6
                read(line(pos+1:), *) mesh%n_groups
            else if (index(line, 'Dims') > 0) then
                pos = index(line, 'Dims') + 4
                read(line(pos+1:), *) mesh%dim
            else if (index(line, '$1D_Patch_Description_Start') > 0) then
                edge_count = edge_count + 1
            else if (index(line, 'POINTS') == 1) then
                read(line(7:), *) mesh%n_nodes
            else if (index(line, 'control_points:') > 0) then
                pos = index(line, ':')
                read(line(pos+1:), *) k
                if (k > max_cp) max_cp = k
            else if (index(line, 'KnotVector') > 0) then
                pos = index(line, ':')
                if (pos == 0) then
                    if (index(line, 'KnotVectorXi') > 0) then
                        pos = index(line, 'KnotVectorXi') + 11
                    else if (index(line, 'KnotVectorEta') > 0) then
                        pos = index(line, 'KnotVectorEta') + 12
                    else
                        pos = index(line, 'KnotVector') + 9
                    end if
                end if
                read(line(pos+1:), *, iostat=iostatus) k
                if (iostatus == 0 .and. k > max_knots) max_knots = k
            end if
        end do

        mesh%n_edges   = edge_count
        mesh%n_faces_per_elem = 4
        p = mesh%order
        q = mesh%order

        allocate(mesh%nodes(mesh%n_nodes, mesh%dim), mesh%weights(mesh%n_nodes))
        allocate(mesh%patch_cp_ids(patch_count, max_cp), mesh%patch_material_ids(patch_count))
        allocate(mesh%patch_n_cp_xi(patch_count), mesh%patch_n_cp_eta(patch_count))
        allocate(mesh%patch_n_knots_xi(patch_count), mesh%patch_n_knots_eta(patch_count))
        allocate(mesh%patch_knot_vectors_xi(patch_count, max_knots), mesh%patch_knot_vectors_eta(patch_count, max_knots))
        allocate(mesh%edge_n_cp_xi(edge_count), mesh%edge_n_knots_xi(edge_count))
        allocate(mesh%edge_knot_vectors_xi(edge_count, max_knots))
        allocate(mesh%edges(edge_count, max_cp), mesh%boundary_ids(edge_count))

        mesh%patch_cp_ids = 0
        mesh%patch_knot_vectors_xi = 0.0
        mesh%patch_knot_vectors_eta = 0.0
        mesh%edges = 0
        mesh%boundary_ids = 0

        rewind(unit)
    
        jj = 0 ! 2D Patch counter
        ii = 0 ! 1D Patch counter
        
        do
            read(unit, '(A)', iostat=iostatus) line
            if (iostatus /= 0) exit
            line = adjustl(line)

            if (len_trim(line) == 0 .or. line(1:1) == '!') cycle

            if (index(line, 'POINTS') == 1) then
                do k = 1, mesh%n_nodes
                    read(unit, *) mesh%nodes(k,1), mesh%nodes(k,2), dummy_coord, mesh%weights(k)
                end do

            else if (index(line, '$1D_Patch_Description_Start') > 0) then
                ii = ii + 1
                call parse_patch_block(unit, mesh%edges(ii,:), mesh%boundary_ids(ii), .false., &
                    n_cp=mesh%edge_n_cp_xi(ii), out_knots_xi=mesh%edge_knot_vectors_xi(ii,:), &
                    n_k_xi=mesh%edge_n_knots_xi(ii))
            

            else if (index(line, '$2D_Patch_Description_Start') > 0) then
                jj = jj + 1
                call parse_patch_block(unit, mesh%patch_cp_ids(jj,:), mesh%patch_material_ids(jj), .true., &
                    out_knots_xi=mesh%patch_knot_vectors_xi(jj,:), out_knots_eta=mesh%patch_knot_vectors_eta(jj,:), &
                    n_k_xi=mesh%patch_n_knots_xi(jj), n_k_eta=mesh%patch_n_knots_eta(jj))
                mesh%patch_n_cp_xi(jj) = mesh%patch_n_knots_xi(jj) - p - 1
                mesh%patch_n_cp_eta(jj) = mesh%patch_n_knots_eta(jj) - q - 1
            end if
        end do
        
        close(unit)

        ! --- Subdivision into Knot Spans (Elements) ---
        n_total_elems = 0
        do p_idx = 1, patch_count
            p = mesh%order
            ! Spans with basis support exist between knot indices p+1 and n_knots - p - 1
            do span_u = p + 1, mesh%patch_n_knots_xi(p_idx) - p - 1
                if (mesh%patch_knot_vectors_xi(p_idx, span_u+1) > mesh%patch_knot_vectors_xi(p_idx, span_u) + dp_EPSILON) then
                    do span_v = p + 1, mesh%patch_n_knots_eta(p_idx) - p - 1
                        if (mesh%patch_knot_vectors_eta(p_idx, span_v+1) > mesh%patch_knot_vectors_eta(p_idx, span_v) + dp_EPSILON) then
                            n_total_elems = n_total_elems + 1
                        end if
                    end do
                end if
            end do
        end do

        mesh%n_elems = n_total_elems
        allocate(mesh%elems(mesh%n_elems, (p+1)*(q+1)))
        allocate(mesh%material_ids(mesh%n_elems))
        allocate(mesh%elem_u_min(mesh%n_elems), mesh%elem_u_max(mesh%n_elems))
        allocate(mesh%elem_v_min(mesh%n_elems), mesh%elem_v_max(mesh%n_elems))
        allocate(mesh%elem_patch_id(mesh%n_elems))
        allocate(mesh%elem_map_to_id(patch_count, max_knots, max_knots)) ! Max_knots is a safe upper bound for span indices
        allocate(mesh%elem_span_indices(2, mesh%n_elems))

        k = 0
        do p_idx = 1, patch_count
            do span_u = p + 1, mesh%patch_n_knots_xi(p_idx) - p - 1
                if (mesh%patch_knot_vectors_xi(p_idx, span_u+1) <= mesh%patch_knot_vectors_xi(p_idx, span_u) + dp_EPSILON) cycle
                do span_v = p + 1, mesh%patch_n_knots_eta(p_idx) - p - 1
                    if (mesh%patch_knot_vectors_eta(p_idx, span_v+1) <= mesh%patch_knot_vectors_eta(p_idx, span_v) + dp_EPSILON) cycle
                    
                    k = k + 1
                    mesh%elem_patch_id(k) = p_idx
                    mesh%elem_span_indices(1, k) = span_u
                    mesh%elem_span_indices(2, k) = span_v
                    mesh%elem_map_to_id(p_idx, span_u, span_v) = k ! Populate lookup table
                    mesh%material_ids(k) = mesh%patch_material_ids(p_idx)
                    mesh%elem_u_min(k)   = mesh%patch_knot_vectors_xi(p_idx, span_u)
                    mesh%elem_u_max(k)   = mesh%patch_knot_vectors_xi(p_idx, span_u+1)
                    mesh%elem_v_min(k)   = mesh%patch_knot_vectors_eta(p_idx, span_v)
                    mesh%elem_v_max(k)   = mesh%patch_knot_vectors_eta(p_idx, span_v+1)

                    cp_idx = 0
                    do jj = span_v - p, span_v
                        do ii = span_u - p, span_u
                            cp_idx = cp_idx + 1
                            mesh%elems(k, cp_idx) = mesh%patch_cp_ids(p_idx, (jj-1)*mesh%patch_n_cp_xi(p_idx) + ii)
                        end do
                    end do
                end do
            end do
        end do
        
        call write_mesh_to_files(mesh)

    end subroutine read_asmg_mesh

    subroutine parse_patch_block(u, out_cp, out_id2, is_2d, n_cp, out_knots_xi, out_knots_eta, n_k_xi, n_k_eta)
        integer, intent(in) :: u
        integer, intent(inout) :: out_cp(:)
        integer, intent(out)   :: out_id2
        logical, intent(in)    :: is_2d
        integer, intent(out), optional :: n_cp
        real(dp), intent(inout), optional :: out_knots_xi(:), out_knots_eta(:)
        integer, intent(out), optional :: n_k_xi, n_k_eta
        
        character(len=1024) :: l
        integer :: ios, n_val, m, pos
        
        ! Initialize defaults
        out_id2 = -1
        
        do
            read(u, '(A)', iostat=ios) l
            if (ios /= 0) exit
            l = adjustl(l)
            
            if (index(l, '$1D_Patch_Description_End') > 0 .or. &
                index(l, '$2D_Patch_Description_End') > 0) exit
            
            if (l == '' .or. l(1:1) == '!') cycle

            if (index(l, 'Material_ID:') > 0) then
                pos = index(l, ':')
                ! Check if ID is on the same line or next line
                if (len_trim(l(pos+1:)) > 0) then
                    read(l(pos+1:), *, iostat=ios) out_id2
                else
                    read(u, *, iostat=ios) out_id2
                end if
            else if (index(l, 'BC:') > 0) then
                pos = index(l, ':')
                ! Check if ID is on the same line or next line
                if (len_trim(l(pos+1:)) > 0) then
                    read(l(pos+1:), *, iostat=ios) out_id2
                else
                    read(u, *, iostat=ios) out_id2
                end if
            else if (index(l, 'control_points:') > 0) then
                pos = index(l, ':') + 1
                read(l(pos:), *, iostat=ios) n_val
                if (present(n_cp)) n_cp = n_val
                if (ios == 0 .and. n_val > 0) then
                    read(u, *, iostat=ios) (out_cp(m), m=1, n_val)
                    if (ios == 0) out_cp(1:n_val) = out_cp(1:n_val) + 1
                end if
            else if (index(l, 'KnotVectorXi') > 0 .and. present(out_knots_xi)) then
                pos = index(l, ':')
                if (pos == 0) pos = index(l, 'KnotVectorXi') + 11
                read(l(pos+1:), *, iostat=ios) n_val
                if (ios == 0 .and. n_val > 0) then
                    read(u, *, iostat=ios) (out_knots_xi(m), m=1, n_val)
                    if (present(n_k_xi)) n_k_xi = n_val
                end if
            else if (index(l, 'KnotVectorEta') > 0 .and. present(out_knots_eta)) then
                pos = index(l, ':')
                if (pos == 0) pos = index(l, 'KnotVectorEta') + 12
                read(l(pos+1:), *, iostat=ios) n_val
                if (ios == 0 .and. n_val > 0) then
                    read(u, *, iostat=ios) (out_knots_eta(m), m=1, n_val)
                    if (present(n_k_eta)) n_k_eta = n_val
                end if
            else if (index(l, 'KnotVector') > 0 .and. present(out_knots_xi)) then
                pos = index(l, ':')
                if (pos == 0) pos = index(l, 'KnotVector') + 9
                read(l(pos+1:), *, iostat=ios) n_val
                if (ios == 0 .and. n_val > 0) then
                    read(u, *, iostat=ios) (out_knots_xi(m), m=1, n_val)
                    if (present(n_k_xi)) n_k_xi = n_val
                    if (is_2d .and. present(out_knots_eta) .and. present(n_k_eta)) then
                         out_knots_eta(1:n_val) = out_knots_xi(1:n_val)
                         n_k_eta = n_val
                    end if
                end if
            end if
        end do
    end subroutine parse_patch_block

    subroutine write_mesh_to_files(mesh)
        type(t_mesh), intent(in) :: mesh
        integer :: u, ii

        ! 1. Export Nodes and Weights
        open(newunit=u, file='../output/nodes.dat', status='replace')
        write(u, '(A)') "# ID | X | Y | Weight"
        do ii = 1, mesh%n_nodes
            write(u, '(I8, 3F15.8)') ii, mesh%nodes(ii,1), mesh%nodes(ii,2), mesh%weights(ii)
        end do
        close(u)

        ! 2. Export Elements (Control Points)
        open(newunit=u, file='../output/elements.dat', status='replace')
        write(u, '(A)') "# Element ID | Control Point IDs..."
        do ii = 1, mesh%n_elems
            write(u, '(I8, " : ", 500I8)') ii, pack(mesh%elems(ii,:), mesh%elems(ii,:) /= 0)
        end do
        close(u)

        ! 3. Export Edges (1D Patches)
        open(newunit=u, file='../output/edges.dat', status='replace')
        write(u, '(A)') "# Edge ID | Boundary ID | Control Point IDs..."
        do ii = 1, mesh%n_edges
            write(u, '(I8, " BC:", I4, " : ", 500I8)') ii, mesh%boundary_ids(ii), &
                  pack(mesh%edges(ii,:), mesh%edges(ii,:) /= 0)
        end do
        close(u)

        ! 4. Export Materials and Pin IDs
        open(newunit=u, file='../output/materials.dat', status='replace')
        write(u, '(A)') "# Elem ID | Material ID"
        do ii = 1, mesh%n_elems
            write(u, '(2I10)') ii, mesh%material_ids(ii)
        end do
        close(u)

        ! 5. Export Parametric Space
        open(newunit=u, file='../output/knotspans.dat', status='replace')
        write(u, '(A)') "# Elem ID | U_range | V_range"
        do ii = 1, mesh%n_elems
            write(u, '(1I10, 4F12.6)') ii, &
                  mesh%elem_u_min(ii), mesh%elem_u_max(ii), &
                  mesh%elem_v_min(ii), mesh%elem_v_max(ii)
        end do
        close(u)

    end subroutine write_mesh_to_files

    subroutine InsertKnot1D(p, knots_old, cp_old, weights_old, u_new, knots_ref, cp_ref, weights_ref)
        integer, intent(in) :: p
        real(dp), intent(in) :: knots_old(:), cp_old(:,:), weights_old(:), u_new
        real(dp), allocatable, intent(out) :: knots_ref(:), cp_ref(:,:), weights_ref(:)
        
        integer :: n, k, i, dim
        real(dp) :: alpha
        
        n = size(cp_old, 2)
        dim = size(cp_old, 1)
        
        ! 1. Find the span where the new knot sits
        k = 0
        do i = 1, size(knots_old)-1
            if (u_new >= knots_old(i) .and. u_new < knots_old(i+1)) then
                k = i
                exit
            end if
        end do
        if (k == 0) return

        allocate(knots_ref(size(knots_old)+1))
        allocate(cp_ref(dim, n+1))
        allocate(weights_ref(n+1))

        ! 2. New knot vector
        knots_ref(1:k) = knots_old(1:k)
        knots_ref(k+1) = u_new
        knots_ref(k+2:) = knots_old(k+1:)

        ! 3. New Control Points and Weights
        do i = 1, n + 1
            if (i <= k - p) then
                alpha = 1.0_dp
            else if (i >= k + 1) then
                alpha = 0.0_dp
            else
                alpha = (u_new - knots_old(i)) / (knots_old(i + p) - knots_old(i))
            end if
            
            if (i == 1) then
                cp_ref(:, i) = cp_old(:, i)
                weights_ref(i) = weights_old(i)
            else if (i == n + 1) then
                cp_ref(:, i) = cp_old(:, i-1)
                weights_ref(i) = weights_old(i-1)
            else
                cp_ref(:, i) = alpha * cp_old(:, i) * weights_old(i) + (1.0_dp - alpha) * cp_old(:, i-1) * weights_old(i-1)
                weights_ref(i) = alpha * weights_old(i) + (1.0_dp - alpha) * weights_old(i-1)
                cp_ref(:, i) = cp_ref(:, i) / weights_ref(i) ! Project back from homogeneous space
            end if
        end do
    end subroutine InsertKnot1D

end module m_asmg
