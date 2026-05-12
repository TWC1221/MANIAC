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
        type(t_mesh), allocatable, intent(inout) :: mesh

        integer :: unit, ios, ii, jj, kk, k, p_idx, s_u, s_v, s_w, cp_idx, n_tot
        character(len=1024) :: line
        integer :: n_p, n_s, max_k, p, q, r, nxi, neta

        if (.not. allocated(mesh)) allocate(mesh)
        open(newunit=unit, file=filepath, status='old', action='read')
        mesh%order = 1; mesh%n_nodes = 0; n_p = 0; n_s = 0; max_k = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (index(line, '$3D_Patch_Description_Start') > 0) n_p = n_p + 1
            if (index(line, '$2D_Patch_Description_Start') > 0) n_s = n_s + 1
            if (index(line, 'POINTS') == 1) read(line(7:), *) mesh%n_nodes
            ii = scan(line, ': ')
            if (ii > 0) then
                if (index(line, 'PolyOrder') > 0) read(line(ii+1:), *) mesh%order
                if (index(line, 'Groups') > 0)    read(line(ii+1:), *) mesh%n_groups
                if (index(line, 'Dims') > 0)      read(line(ii+1:), *) mesh%dim
                if (index(line, 'KnotVector') > 0) then
                    read(line(ii+1:), *, iostat=ios) jj
                    if (ios == 0 .and. jj > max_k) max_k = jj
                end if
            end if
        end do
        
        ! Safety: ensure max_k is sufficient for the knot vector indexing
        max_k = max_k + 5 

        p = mesh%order; q = p; r = p
        mesh%n_faces_per_elem = 6
        allocate(mesh%nodes(mesh%n_nodes, 3), mesh%weights(mesh%n_nodes))
        allocate(mesh%patches(n_p), mesh%surfaces(n_s))
        do k = 1, n_p; mesh%patches(k)%face_to_surface = 0; end do
        rewind(unit)

        ii = 0; jj = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '!') cycle

            if (index(line, 'POINTS') == 1) then
                do k = 1, mesh%n_nodes
                    read(unit, *) mesh%nodes(k, 1:3), mesh%weights(k)
                end do
            else if (index(line, '$2D_Patch_Description_Start') > 0) then
                ii = ii + 1; call parse_block(unit, mesh%surfaces(ii)%cp_ids, mesh%surfaces(ii)%bc_id, &
                    mesh%surfaces(ii)%knots_xi, mesh%surfaces(ii)%knots_eta)
            else if (index(line, '$3D_Patch_Description_Start') > 0) then
                jj = jj + 1
                call parse_block(unit, mesh%patches(jj)%cp_ids, mesh%patches(jj)%material_id, &
                    mesh%patches(jj)%knots_xi, mesh%patches(jj)%knots_eta, mesh%patches(jj)%knots_zeta)
            end if
        end do
        close(unit)

        ! --- Subdivision into Knot Spans (Elements) ---
        n_tot = 0
        do p_idx = 1, n_p
            associate(ptch => mesh%patches(p_idx))
            do s_u = p + 1, size(ptch%knots_xi) - p - 1
                if (ptch%knots_xi(s_u+1) > ptch%knots_xi(s_u) + dp_EPSILON) then
                    do s_v = q + 1, size(ptch%knots_eta) - q - 1
                        if (ptch%knots_eta(s_v+1) > ptch%knots_eta(s_v) + dp_EPSILON) then
                            do s_w = r + 1, size(ptch%knots_zeta) - r - 1
                                if (ptch%knots_zeta(s_w+1) > ptch%knots_zeta(s_w) + dp_EPSILON) n_tot = n_tot + 1
                            end do
                        end if
                    end do
                end if
            end do
            end associate
        end do

        mesh%n_elems = n_tot; k = (p+1)**3
        allocate(mesh%elems(n_tot, k), mesh%material_ids(n_tot), mesh%elem_u_min(n_tot), &
                 mesh%elem_u_max(n_tot), mesh%elem_v_min(n_tot), mesh%elem_v_max(n_tot), &
                 mesh%elem_w_min(n_tot), mesh%elem_w_max(n_tot), mesh%elem_patch_id(n_tot), &
                 mesh%elem_span_indices(3, n_tot), &
                 mesh%elem_map_to_id(n_p, max_k, max_k, max_k))

        k = 0
        do p_idx = 1, n_p
            associate(ptch => mesh%patches(p_idx))
            nxi = size(ptch%knots_xi) - p - 1; neta = size(ptch%knots_eta) - q - 1
            do s_u = p + 1, size(ptch%knots_xi) - p - 1
                if (ptch%knots_xi(s_u+1) <= ptch%knots_xi(s_u) + dp_EPSILON) cycle
                do s_v = q + 1, size(ptch%knots_eta) - q - 1
                    if (ptch%knots_eta(s_v+1) <= ptch%knots_eta(s_v) + dp_EPSILON) cycle
                    do s_w = r + 1, size(ptch%knots_zeta) - r - 1
                        if (ptch%knots_zeta(s_w+1) <= ptch%knots_zeta(s_w) + dp_EPSILON) cycle

                        k = k + 1
                        mesh%elem_patch_id(k) = p_idx; mesh%elem_span_indices(1, k) = s_u
                        mesh%elem_span_indices(2, k) = s_v; mesh%elem_span_indices(3, k) = s_w
                        mesh%elem_map_to_id(p_idx, s_u, s_v, s_w) = k
                        
                        mesh%material_ids(k) = ptch%material_id
                        mesh%elem_u_min(k) = ptch%knots_xi(s_u); mesh%elem_u_max(k) = ptch%knots_xi(s_u+1)
                        mesh%elem_v_min(k) = ptch%knots_eta(s_v); mesh%elem_v_max(k) = ptch%knots_eta(s_v+1)
                        mesh%elem_w_min(k) = ptch%knots_zeta(s_w); mesh%elem_w_max(k) = ptch%knots_zeta(s_w+1)

                        cp_idx = 0
                        do kk = s_w - r, s_w
                            do jj = s_v - q, s_v
                                do ii = s_u - p, s_u
                                    cp_idx = cp_idx + 1
                                    mesh%elems(k, cp_idx) = ptch%cp_ids((kk-1)*nxi*neta + (jj-1)*nxi + ii)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            end associate
        end do
        call write_mesh_to_files(mesh)
    end subroutine read_asmg_mesh

    subroutine parse_block(u, cp, id, k1, k2, k3)
        integer, intent(in) :: u
        integer, allocatable, intent(out) :: cp(:)
        integer, intent(out) :: id
        real(dp), allocatable, intent(out) :: k1(:), k2(:)
        real(dp), allocatable, intent(out), optional :: k3(:)
        
        character(len=1024) :: l, key
        integer :: n, ios
        real(dp), allocatable :: temp_k(:)

        id = -1
        do
            read(u, '(A)', iostat=ios) l
            if (ios /= 0 .or. index(l, '_End') > 0) exit
            l = adjustl(l)
            if (l == '' .or. l(1:1) == '!') cycle

            ! Use internal read to split keyword from immediate value (if any)
            read(l, *, iostat=ios) key, n 

            if (index(l, 'control_points') > 0) then
                if (ios /= 0) read(u, *) n  ! If value wasn't on same line, read next
                allocate(cp(n)); read(u, *) cp; cp = cp + 1

            else if (index(l, 'KnotVector_Xi') > 0) then
                if (ios /= 0) read(u, *) n
                if (allocated(k1)) deallocate(k1)
                allocate(k1(n)); read(u, *) k1
            else if (index(l, 'KnotVector_Eta') > 0) then
                if (ios /= 0) read(u, *) n
                if (allocated(k2)) deallocate(k2)
                allocate(k2(n)); read(u, *) k2
            else if (index(l, 'KnotVector_Zeta') > 0) then
                if (ios /= 0) read(u, *) n
                if (present(k3)) then
                    if (allocated(k3)) deallocate(k3)
                    allocate(k3(n)); read(u, *) k3
                end if
            else if (index(l, 'KnotVector') > 0) then
                if (ios /= 0) read(u, *) n
                allocate(temp_k(n)); read(u, *) temp_k
                if (.not. allocated(k1)) allocate(k1(n), source=temp_k)
                if (.not. allocated(k2)) allocate(k2(n), source=temp_k)
                if (present(k3)) then; if (.not. allocated(k3)) allocate(k3(n), source=temp_k); end if
                deallocate(temp_k)

            else if (index(l, 'Material_ID') > 0 .or. index(l, 'BC') > 0) then
                if (ios == 0) then; id = n; else; read(u, *) id; end if
            end if
        end do
    end subroutine parse_block

    subroutine write_mesh_to_files(mesh)
        type(t_mesh), intent(in) :: mesh
        integer :: u, ii

        open(newunit=u, file='../output/nodes.dat', status='replace')
        write(u, '(A)') "# ID | X | Y | Z | Weight"
        do ii = 1, mesh%n_nodes; write(u, '(I8, 4F15.8)') ii, mesh%nodes(ii,1:3), mesh%weights(ii); end do
        close(u)

        open(newunit=u, file='../output/elements.dat', status='replace')
        write(u, '(A)') "# Element ID | Control Point IDs..."
        do ii = 1, mesh%n_elems; write(u, '(I8, " : ", 500I8)') ii, mesh%elems(ii, :); end do
        close(u)

        open(newunit=u, file='../output/edges.dat', status='replace')
        write(u, '(A)') "# Edge ID | Boundary ID | Control Point IDs..."
        do ii = 1, size(mesh%surfaces); write(u, '(I8, " BC:", I4, " : ", 500I8)') ii, mesh%surfaces(ii)%bc_id, mesh%surfaces(ii)%cp_ids; end do
        close(u)

        open(newunit=u, file='../output/materials.dat', status='replace')
        write(u, '(A)') "# Elem ID | Material ID"
        do ii = 1, mesh%n_elems; write(u, '(2I10)') ii, mesh%material_ids(ii); end do
        close(u)

        open(newunit=u, file='../output/knotspans.dat', status='replace')
        write(u, '(A)') "# Elem ID | U_range | V_range | W_range"
        do ii = 1, mesh%n_elems
            write(u, '(1I10, 6F12.6)') ii, &
                  mesh%elem_u_min(ii), mesh%elem_u_max(ii), &
                  mesh%elem_v_min(ii), mesh%elem_v_max(ii), &
                  mesh%elem_w_min(ii), mesh%elem_w_max(ii)
        end do
        close(u)

    end subroutine write_mesh_to_files

end module m_asmg
