module m_GMSH
    use m_constants
    use m_types
    use m_quadrature
    implicit none
    private
    public :: ParseMesh, ParsePinPowers, derive_case_nametag

contains

    pure function derive_case_nametag(filename) result(tag)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: tag
        integer :: pos, dot_pos

        ! 1. Strip the path (keep only the filename)
        pos = index(filename, '/', back=.true.)
        if (pos == 0) pos = index(filename, '\', back=.true.)

        if (pos > 0) then
            tag = filename(pos+1:)
        else
            tag = filename
        end if

        ! 2. Strip the .vtk extension if it exists
        dot_pos = index(tag, '.vtk', back=.true.)
        if (dot_pos > 0) then
            tag = tag(:dot_pos-1)
        end if
    end function derive_case_nametag

    subroutine ParseMesh(filename, FE, mesh, is_SEM)
        character(len=*), intent(in)    :: filename
        logical,        intent(in)      :: is_SEM
        type(t_finite), intent(inout)   :: FE
        type(t_mesh),   intent(out)     :: mesh

        integer :: unit, ios
        integer :: nnode, nelem, total_ints, target_p
        integer :: i, k, npe, start, max_npe
        character(len=256) :: line, key, dtype, output_dir

        real(dp), allocatable :: nodes3(:,:)
        integer,  allocatable :: connectivity(:)
        integer,  allocatable :: elem_ptr(:)
        integer,  allocatable :: elem_type(:)
        integer,  allocatable :: material_all(:)
        integer,  allocatable :: pin_all(:)

        integer :: n_hexs, n_quads, n_edge_nodes, n_vol_nodes
        integer :: ie, ib

        call timer_start('Total Execution')
        call timer_start('Initialisation')
        call timer_start('I/O  - Parse Mesh')

        output_dir = "../output/" // derive_case_nametag(filename)
        call execute_command_line("mkdir -p '" // trim(output_dir) // "/mesh_cache'")

        print*, '>>> Input Mesh: '//trim(filename)

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'm_GMSH: cannot open file'

        ! --- POINTS section ---
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) stop 'm_GMSH: POINTS section not found'
            if (index(line, 'POINTS') > 0) exit
        end do
        
        read(line,*) key, nnode, dtype
        allocate(nodes3(nnode,3))
        do i = 1, nnode
            read(unit,*) nodes3(i,:)
        end do

        ! --- CELLS section ---
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) stop 'm_GMSH: CELLS section not found'
            if (index(line, 'CELLS') > 0) exit
        end do
        read(line,*) key, nelem, total_ints
        allocate(connectivity(total_ints))
        allocate(elem_ptr(nelem+1))
        read(unit,*) connectivity

        k = 1
        elem_ptr(1) = 1
        max_npe = 0
        do i = 1, nelem
            npe = connectivity(k)
            max_npe = max(max_npe, npe)
            elem_ptr(i+1) = elem_ptr(i) + 1 + npe
            k = elem_ptr(i+1)
        end do

        ! --- CELL_TYPES section ---
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) stop 'm_GMSH: CELL_TYPES section not found'
            if (index(line, 'CELL_TYPES') > 0) exit
        end do
        allocate(elem_type(nelem))
        do i = 1, nelem
            read(unit,*) elem_type(i)
        end do

        ! --- CELL_DATA (Materials) ---
        allocate(material_all(nelem)); material_all = 1
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'CELL_DATA') > 0) exit
        end do

        if (ios == 0) then
            read(unit,*) ! SCALARS...
            read(unit,*) ! LOOKUP_TABLE...
            do i = 1, nelem
                read(unit,*,iostat=ios) material_all(i)
                if (ios /= 0) material_all(i) = 1
            end do
        end if

        ! --- PinID Data ---
        allocate(pin_all(nelem)); pin_all = 0
        do 
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'SCALARS PinID') > 0) exit
        end do

        if (ios == 0) then
            read(unit,*) ! LOOKUP_TABLE...
            do i = 1, nelem
                read(unit,*,iostat=ios) pin_all(i)
                if (ios /= 0) pin_all(i) = 0
            end do
        end if
        close(unit)

        target_p = FE%order
        mesh%n_nodes = nnode
        allocate(mesh%nodes(mesh%n_nodes, 3)) 

        mesh%nodes = nodes3
        if (allocated(nodes3)) deallocate(nodes3)

        ! Check dimensions
        mesh%dim = 3
        n_hexs  = count(elem_type == 12 .or. elem_type == 29) ! Hex or TriQuadHex
        n_quads = count(elem_type == 9  .or. elem_type == 28) ! Quad or BiQuad

        mesh%n_elems = n_hexs
        mesh%n_edges = n_quads ! In 3D, 'edges' holds boundary faces

        n_vol_nodes  = 0
        n_edge_nodes = 0 ! Boundary nodes

        do i = 1, nelem
            npe = connectivity(elem_ptr(i))
            if (elem_type(i) == 12 .or. elem_type(i) == 29) n_vol_nodes  = max(n_vol_nodes, npe)
            if (elem_type(i) == 9  .or. elem_type(i) == 28) n_edge_nodes = max(n_edge_nodes, npe)
        end do

        FE%n_basis = n_vol_nodes
        FE%order = nint(real(FE%n_basis, dp)**(1.0_dp/3.0_dp)) - 1
        FE%n_nodes_per_face = (FE%order + 1)**2

        allocate(mesh%elems(mesh%n_elems, n_vol_nodes))
        allocate(mesh%edges(mesh%n_edges, n_edge_nodes))
        allocate(mesh%material_ids(mesh%n_elems))
        allocate(mesh%pin_ids(mesh%n_elems))
        allocate(mesh%boundary_ids(mesh%n_edges))

        mesh%edges = 0
        mesh%boundary_ids = 0

        ie = 0; ib = 0
        do i = 1, nelem
            start = elem_ptr(i)
            npe   = connectivity(start)

            select case (elem_type(i))
            case (12, 29) ! VTK_HEXAHEDRON
                ie = ie + 1
                mesh%material_ids(ie) = material_all(i)
                mesh%pin_ids(ie) = pin_all(i)
                ! Map VTK to FE Tensor (x-fastest, y, z)
                ! VTK: 0(BLB), 1(BRB), 2(TRB), 3(TLB), 4(BLT), 5(BRT), 6(TRT), 7(TLT)
                ! FE:  1(BLB), 2(BRB), 3(TLB), 4(TRB), 5(BLT), 6(BRT), 7(TLT), 8(TRT)
                if (npe == 8) then
                    mesh%elems(ie, 1) = connectivity(start + 1) + 1
                    mesh%elems(ie, 2) = connectivity(start + 2) + 1
                    mesh%elems(ie, 3) = connectivity(start + 4) + 1 ! Swap 3 & 4
                    mesh%elems(ie, 4) = connectivity(start + 3) + 1 
                    mesh%elems(ie, 5) = connectivity(start + 5) + 1
                    mesh%elems(ie, 6) = connectivity(start + 6) + 1
                    mesh%elems(ie, 7) = connectivity(start + 8) + 1 ! Swap 7 & 8
                    mesh%elems(ie, 8) = connectivity(start + 7) + 1
                else
                    mesh%elems(ie, :) = connectivity(start+1:start+npe) + 1
                end if

            case (9, 28) ! VTK_QUAD, VTK_BIQUADRATIC_QUAD
                ! 3D Boundary Face
                ib = ib + 1
                mesh%boundary_ids(ib) = material_all(i)
                if (npe == 4) then
                    mesh%edges(ib, 1) = connectivity(start + 1) + 1
                    mesh%edges(ib, 2) = connectivity(start + 2) + 1
                    mesh%edges(ib, 4) = connectivity(start + 3) + 1 ! Swap for tensor (TR <-> TL)
                    mesh%edges(ib, 3) = connectivity(start + 4) + 1 
                else
                    mesh%edges(ib, :) = connectivity(start+1:start+npe) + 1
                end if
            end select
        end do

        if (target_p > 1 .and. FE%order == 1) then
            write(*,*) ">>> Upgrading Trilinear Hexahedral Mesh to Order ", trim(int_to_str(target_p))
            call ConvertMeshToHighOrder(mesh, target_p, FE, is_SEM)
        end if

        call write_nodes(mesh%nodes, output_dir)
        call write_elements(nelem, connectivity, elem_ptr, elem_type, output_dir)
        call write_edges(mesh%n_edges, mesh%edges, mesh%boundary_ids, output_dir)
        call write_materials(mesh%n_elems, mesh%material_ids, output_dir)
        call write_pins(mesh%n_elems, mesh%pin_ids, output_dir)

        if (allocated(connectivity))  deallocate(connectivity)
        if (allocated(elem_ptr))      deallocate(elem_ptr)
        if (allocated(elem_type))     deallocate(elem_type)
        if (allocated(material_all))  deallocate(material_all)
        if (allocated(pin_all))       deallocate(pin_all)

        call timer_stop('I/O  - Parse Mesh')
    end subroutine ParseMesh

    subroutine ConvertMeshToHighOrder(mesh, p, FE, is_SEM)
        type(t_mesh),   intent(inout) :: mesh
        integer,        intent(in)    :: p
        type(t_finite), intent(inout) :: FE
        logical,        intent(in)    :: is_SEM
        
        integer :: n_lin_nodes, n_elems, n_edges_unique
        integer :: i, j, k, n1, n2
        integer :: n_new_nodes_per_edge, n_new_nodes_per_face, n_new_nodes_per_elem_internal
        integer :: node_count, idx
        integer, allocatable :: unique_edge_list(:,:)
        type(t_edge_sort), allocatable :: sort_list(:)
        real(dp), allocatable :: new_nodes(:,:)
        integer, allocatable  :: new_elems(:,:)
        integer :: n_basis_ho, ix, iy, iz
        type(t_quadrature) :: QuadLob
        real(dp) :: tmp_pos(3)

        n_lin_nodes = mesh%n_nodes
        n_elems     = mesh%n_elems

        ! 1. Identify and Sort Edges to find unique ones
        allocate(sort_list(n_elems * 12))
        k = 0
        do i = 1, n_elems
            ! Standard CCW/Lexicographic edge extraction
            call add_edge_to_list(sort_list, k, mesh%elems(i,1), mesh%elems(i,2))
            call add_edge_to_list(sort_list, k, mesh%elems(i,2), mesh%elems(i,4))
            call add_edge_to_list(sort_list, k, mesh%elems(i,4), mesh%elems(i,3))
            call add_edge_to_list(sort_list, k, mesh%elems(i,3), mesh%elems(i,1))
            call add_edge_to_list(sort_list, k, mesh%elems(i,1), mesh%elems(i,5))
            call add_edge_to_list(sort_list, k, mesh%elems(i,2), mesh%elems(i,6))
            call add_edge_to_list(sort_list, k, mesh%elems(i,4), mesh%elems(i,8))
            call add_edge_to_list(sort_list, k, mesh%elems(i,3), mesh%elems(i,7))
            call add_edge_to_list(sort_list, k, mesh%elems(i,5), mesh%elems(i,6))
            call add_edge_to_list(sort_list, k, mesh%elems(i,6), mesh%elems(i,8))
            call add_edge_to_list(sort_list, k, mesh%elems(i,8), mesh%elems(i,7))
            call add_edge_to_list(sort_list, k, mesh%elems(i,7), mesh%elems(i,5))
        end do

        call sort_edges(sort_list, k)
        
        ! Count unique edges
        n_edges_unique = 0
        if (k > 0) n_edges_unique = 1
        do i = 2, k
            if (sort_list(i)%n1 /= sort_list(i-1)%n1 .or. sort_list(i)%n2 /= sort_list(i-1)%n2) then
                n_edges_unique = n_edges_unique + 1
            end if
        end do


        ! Compress to unique list
        allocate(unique_edge_list(n_edges_unique, 2))
        j = 0
        if (k > 0) then
            j = 1
            unique_edge_list(1,1) = sort_list(1)%n1
            unique_edge_list(1,2) = sort_list(1)%n2
        end if
        do i = 2, k
            if (sort_list(i)%n1 /= sort_list(i-1)%n1 .or. sort_list(i)%n2 /= sort_list(i-1)%n2) then
                j = j + 1
                unique_edge_list(j,1) = sort_list(i)%n1
                unique_edge_list(j,2) = sort_list(i)%n2
            end if
        end do

        ! 2. Allocate New Arrays (Using Contiguous 2-by-N Layout)
        n_new_nodes_per_edge = p - 1
        n_new_nodes_per_face = (p - 1)**2
        n_new_nodes_per_elem_internal = (p - 1)**3
        n_basis_ho = (p + 1)**3

        ! Note: For a strictly continuous mesh, face nodes should also be searched and shared.
        mesh%n_nodes = n_lin_nodes + n_edges_unique * n_new_nodes_per_edge + n_elems * (6*n_new_nodes_per_face + n_new_nodes_per_elem_internal)
        allocate(new_nodes(mesh%n_nodes, 3))
        allocate(new_elems(n_elems, n_basis_ho))

        new_nodes = 0.0_dp
        new_elems = 0

        call SpectralLinearQuadrature(QuadLob, p)

        ! Copy linear nodes
        new_nodes(1:n_lin_nodes, :) = mesh%nodes(1:n_lin_nodes, :)

        ! Generate Edge Nodes
        node_count = n_lin_nodes
        do i = 1, n_edges_unique
            n1 = unique_edge_list(i, 1)
            n2 = unique_edge_list(i, 2)
            do j = 1, n_new_nodes_per_edge
                node_count = node_count + 1
                ! Passing (node_count, :)
                if (is_SEM) then
                    call real_interp(new_nodes(node_count, :), mesh%nodes(n1, :), mesh%nodes(n2, :), (QuadLob%xi(j+1) + 1.0_dp)/2.0_dp)
                else
                    call real_interp(new_nodes(node_count, :), mesh%nodes(n1, :), mesh%nodes(n2, :), real(j,dp)/real(p,dp))
                end if
            end do
        end do

        ! 3. Build Connectivity and Internal Nodes
        do i = 1, n_elems
            ! Corners (Lexicographic: BL, BR, TL, TR)
            new_elems(i, 1)                       = mesh%elems(i, 1)
            new_elems(i, p+1)                     = mesh%elems(i, 2)
            new_elems(i, (p+1)*p + 1)             = mesh%elems(i, 3)
            new_elems(i, (p+1)**2)                = mesh%elems(i, 4)
            new_elems(i, 1 + p*(p+1)**2)          = mesh%elems(i, 5)
            new_elems(i, p+1 + p*(p+1)**2)        = mesh%elems(i, 6)
            new_elems(i, (p+1)*p + 1 + p*(p+1)**2)= mesh%elems(i, 7)
            new_elems(i, (p+1)**3)                = mesh%elems(i, 8)

            ! Edge 1: Bottom (BL->BR)
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,1), mesh%elems(i,2), 1, 1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)
            
            ! Edge 2: Right (BR->TR) -> Stride p+1
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,2), mesh%elems(i,4), p+1, p+1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)

            ! Edge 3: Top (TL->TR) -> Stride 1.
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,3), mesh%elems(i,4), (p+1)*p+1, 1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)

            ! Edge 4: Left (BL->TL) -> Stride p+1
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,1), mesh%elems(i,3), 1, p+1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)

            ! Edge 5-8: Mid-height vertical edges
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,1), mesh%elems(i,5), 1, (p+1)**2, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,2), mesh%elems(i,6), p+1, (p+1)**2, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,3), mesh%elems(i,7), (p+1)*p+1, (p+1)**2, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,4), mesh%elems(i,8), (p+1)**2, (p+1)**2, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)

            ! Edge 9-12: Top face edges
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,5), mesh%elems(i,6), 1+p*(p+1)**2, 1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,6), mesh%elems(i,8), p+1+p*(p+1)**2, p+1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,7), mesh%elems(i,8), (p+1)*p+1+p*(p+1)**2, 1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)
            call process_edge_nodes(new_elems(i,:), mesh%elems(i,5), mesh%elems(i,7), 1+p*(p+1)**2, p+1, &
                                    unique_edge_list, n_edges_unique, n_lin_nodes, n_new_nodes_per_edge)

            ! Face Nodes (Simplification: Generate per element)
            do ix = 0, p; do iy = 0, p; do iz = 0, p
                ! Skip corners and edges (already handled)
                k = count([ix==0 .or. ix==p, iy==0 .or. iy==p, iz==0 .or. iz==p])
                if (k /= 1) cycle ! Process only face-interior nodes (where exactly one coord is on boundary)
                
                idx = 1 + ix + iy*(p+1) + iz*(p+1)**2
                if (new_elems(i, idx) == 0) then
                    node_count = node_count + 1
                    new_elems(i, idx) = node_count
                    call get_ho_node_pos(tmp_pos, mesh, i, ix, iy, iz, p, is_SEM, QuadLob)
                    new_nodes(node_count, :) = tmp_pos
                end if
            end do; end do; end do

            ! Internal Nodes
            do iz = 1, p-1
                do iy = 1, p-1
                    do ix = 1, p-1
                        idx = 1 + ix + iy*(p+1) + iz*(p+1)**2
                        node_count = node_count + 1
                        new_elems(i, idx) = node_count
                        call get_ho_node_pos(tmp_pos, mesh, i, ix, iy, iz, p, is_SEM, QuadLob)
                        new_nodes(node_count, :) = tmp_pos
                    end do
                end do
            end do
        end do

        ! Swap pointers and update Finite Element type
        if (allocated(mesh%nodes)) deallocate(mesh%nodes)
        if (allocated(mesh%elems)) deallocate(mesh%elems)
        
        allocate(mesh%nodes(mesh%n_nodes, 3))
        mesh%nodes = new_nodes
        
        allocate(mesh%elems(n_elems, n_basis_ho))
        mesh%elems = new_elems

        FE%order = p
        FE%n_basis = n_basis_ho
        FE%n_nodes_per_face = (p + 1)**2

        if (allocated(QuadLob%xi)) deallocate(QuadLob%xi)
        if (allocated(QuadLob%weights))  deallocate(QuadLob%weights)
        deallocate(new_nodes, new_elems, sort_list, unique_edge_list)
    end subroutine ConvertMeshToHighOrder

    subroutine get_ho_node_pos(pos, mesh, i_el, ix, iy, iz, p, is_SEM, QuadLob)
        real(dp), intent(out) :: pos(3)
        type(t_mesh), intent(in) :: mesh
        integer, intent(in) :: i_el, ix, iy, iz, p
        logical, intent(in) :: is_SEM
        type(t_quadrature), intent(in) :: QuadLob
        real(dp) :: xi, eta, zeta, N_val(8)
        integer :: j

        if (is_SEM) then
            xi  = QuadLob%xi(ix+1)
            eta = QuadLob%xi(iy+1)
            zeta= QuadLob%xi(iz+1)
        else
            xi  = -1.0_dp + 2.0_dp * real(ix, dp) / real(p, dp)
            eta = -1.0_dp + 2.0_dp * real(iy, dp) / real(p, dp)
            zeta= -1.0_dp + 2.0_dp * real(iz, dp) / real(p, dp)
        end if

        ! Trilinear map from Hex-8 corners
        N_val(1) = 0.125_dp * (1.0_dp - xi) * (1.0_dp - eta) * (1.0_dp - zeta)
        N_val(2) = 0.125_dp * (1.0_dp + xi) * (1.0_dp - eta) * (1.0_dp - zeta)
        N_val(3) = 0.125_dp * (1.0_dp - xi) * (1.0_dp + eta) * (1.0_dp - zeta)
        N_val(4) = 0.125_dp * (1.0_dp + xi) * (1.0_dp + eta) * (1.0_dp - zeta)
        N_val(5) = 0.125_dp * (1.0_dp - xi) * (1.0_dp - eta) * (1.0_dp + zeta)
        N_val(6) = 0.125_dp * (1.0_dp + xi) * (1.0_dp - eta) * (1.0_dp + zeta)
        N_val(7) = 0.125_dp * (1.0_dp - xi) * (1.0_dp + eta) * (1.0_dp + zeta)
        N_val(8) = 0.125_dp * (1.0_dp + xi) * (1.0_dp + eta) * (1.0_dp + zeta)

        pos = 0.0_dp
        do j = 1, 8
            pos = pos + N_val(j) * mesh%nodes(mesh%elems(i_el, j), :)
        end do
    end subroutine get_ho_node_pos

    subroutine process_edge_nodes(elem_conn, nA, nB, start_pos, stride, ulist, n_edges, base_n, n_per_edge)
        integer, intent(inout) :: elem_conn(:)
        integer, intent(in) :: nA, nB, start_pos, stride, n_edges, base_n, n_per_edge
        integer, intent(in) :: ulist(:,:)
        
        integer :: mn, mx, idx, j, node_idx, pos
        logical :: reversed

        mn = min(nA, nB)
        mx = max(nA, nB)
        reversed = (nA > nB)

        idx = binary_search_edge(ulist, n_edges, mn, mx)
        
        if (idx > 0) then
             node_idx = base_n + (idx - 1) * n_per_edge
             do j = 1, n_per_edge
                 pos = start_pos + j * stride
                 if (reversed) then
                     elem_conn(pos) = node_idx + (n_per_edge - j + 1)
                 else
                     elem_conn(pos) = node_idx + j
                 end if
             end do
        endif
    end subroutine process_edge_nodes

    function binary_search_edge(list, n, n1, n2) result(idx)
        integer, intent(in) :: list(:,:), n, n1, n2
        integer :: idx, low, high, mid
        low = 1; high = n; idx = -1
        do while (low <= high)
            mid = (low + high) / 2
            if (list(mid,1) < n1 .or. (list(mid,1) == n1 .and. list(mid,2) < n2)) then
                low = mid + 1
            else if (list(mid,1) > n1 .or. (list(mid,1) == n1 .and. list(mid,2) > n2)) then
                high = mid - 1
            else
                idx = mid; exit
            end if
        end do
    end function binary_search_edge

    subroutine add_edge_to_list(list, k, n1, n2)
        type(t_edge_sort), intent(inout) :: list(:)
        integer, intent(inout) :: k
        integer, intent(in) :: n1, n2
        k = k + 1
        list(k)%n1 = min(n1, n2)
        list(k)%n2 = max(n1, n2)
    end subroutine add_edge_to_list

    subroutine sort_edges(list, n)
        type(t_edge_sort), intent(inout) :: list(:)
        integer, intent(in) :: n
        call quicksort(list, 1, n)
    end subroutine sort_edges

    recursive subroutine quicksort(a, first, last)
        type(t_edge_sort), intent(inout) :: a(:)
        integer, intent(in) :: first, last
        integer :: i, j, pivot_n1, pivot_n2
        type(t_edge_sort) :: temp

        if (first < last) then
            pivot_n1 = a(first)%n1
            pivot_n2 = a(first)%n2
            i = first
            j = last
            do while (i < j)
                do while (a(i)%n1 < pivot_n1 .or. (a(i)%n1 == pivot_n1 .and. a(i)%n2 <= pivot_n2))
                    if (i >= last) exit
                    i = i + 1
                end do
                do while (a(j)%n1 > pivot_n1 .or. (a(j)%n1 == pivot_n1 .and. a(j)%n2 > pivot_n2))
                    j = j - 1
                end do
                if (i < j) then
                    temp = a(i); a(i) = a(j); a(j) = temp
                end if
            end do
            temp = a(first); a(first) = a(j); a(j) = temp
            call quicksort(a, first, j-1)
            call quicksort(a, j+1, last)
        end if
    end subroutine quicksort

    subroutine real_interp(res, n1, n2, t)
        real(dp), intent(out) :: res(:)
        real(dp), intent(in)  :: n1(:), n2(:)
        real(dp), intent(in)  :: t
        res = (1.0_dp - t)*n1 + t*n2
    end subroutine real_interp

    subroutine ParsePinPowers(reference_powers, filename)
        real(dp), allocatable, intent(out) :: reference_powers(:)
        character(len=*), intent(in) :: filename

        real(dp), allocatable :: raw(:,:)
        integer :: i, j, unit, ios, n = 34
        integer :: max_pin_id, pin_id, total_nx, total_ny

        allocate(raw(n, n))
        raw = 0.0_dp
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) "Error: Cannot open pin power table: ", trim(filename)
            if (allocated(raw)) deallocate(raw)
            return
        end if

        ! Read upper-triangular matrix row-wise from PINS.txt
        do j = 1, n
            read(unit, *) raw(j, j:n)
        end do
        close(unit)

        ! Reflect upper triangle to lower triangle to make a full matrix
        do j = 1, n
            do i = j + 1, n
                raw(i, j) = raw(j, i)
            end do
        end do

        ! Convert 2D spatial map to a 1D array indexed by pin_id.
        ! assumes a structured pin layout where pin_id = (row_from_bottom-1)*cols + col
        total_nx = n
        total_ny = n
        max_pin_id = total_nx * total_ny
        allocate(reference_powers(0:max_pin_id))
        reference_powers = 0.0_dp

        do j = 1, total_ny ! y-coordinate, or row (from top, 1 to n)
            do i = 1, total_nx ! x-coordinate, or col (from left, 1 to n)
                pin_id = (total_ny - j) * total_nx + i
                if (pin_id <= max_pin_id) then
                    reference_powers(pin_id) = raw(j, i)
                end if
            end do
        end do

        deallocate(raw)
    end subroutine ParsePinPowers

    subroutine write_edges(n_lines, edges, edge_mats, out_dir)
        integer, intent(in) :: n_lines, edges(:,:), edge_mats(:)
        character(len=*), intent(in) :: out_dir
        integer :: i, jj, unit
        open(newunit=unit, file=trim(out_dir)//"/mesh_cache/edges.dat", status="replace")
        write(unit,'(A)') '# Edge_ID    BC_ID      Nodes...'
        do i = 1, n_lines
            write(unit,'(I8,I12)',advance="no") i, edge_mats(i)
            do jj = 1, size(edges, 2)
                if (edges(i,jj) > 0) write(unit,'(I8)',advance="no") edges(i,jj)
            end do
            write(unit,*)
        end do
        close(unit)
    end subroutine write_edges

    subroutine write_nodes(nodes3, out_dir)
        real(dp), intent(in) :: nodes3(:,:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/mesh_cache/nodes.dat", status="replace")
        write(unit,'(A)') '# id          x                      y                      z'
        do i = 1, size(nodes3, 1)
            write(unit,'(I8,3ES25.15)') i, nodes3(i,:)
        end do
        close(unit)
    end subroutine write_nodes

    subroutine write_elements(nelem, connectivity, elem_ptr, elem_type, out_dir)
        integer, intent(in) :: nelem, connectivity(:), elem_ptr(:), elem_type(:)
        character(len=*), intent(in) :: out_dir
        integer :: i, jj, npe, start, unit
        open(newunit=unit, file=trim(out_dir)//"/mesh_cache/elements.dat", status="replace")
        write(unit,'(A)') '# id     vtk_type        nodes...'
        do i = 1, nelem
            start = elem_ptr(i)
            npe   = connectivity(start)
            write(unit,'(I8,I12)',advance="no") i, elem_type(i)
            do jj = 1, npe
                write(unit,'(I8)',advance="no") connectivity(start+jj) + 1
            end do
            write(unit,*)
        end do
        close(unit)
    end subroutine write_elements

    subroutine write_materials(n_elems, mats, out_dir)
        integer, intent(in) :: n_elems, mats(:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/mesh_cache/mats.dat", status="replace")
        write(unit,'(A)') '# Elem_ID    Mat_ID'
        do i = 1, n_elems
            write(unit,'(I8, I10)') i, mats(i)
        end do
        close(unit)
    end subroutine write_materials

    subroutine write_pins(n_elems, pins, out_dir)
        integer, intent(in) :: n_elems, pins(:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/mesh_cache/pins.dat", status="replace")
        write(unit,'(A)') '# Elem_ID    Pin_ID'
        do i = 1, n_elems
            write(unit,'(I8, I10)') i, pins(i)
        end do
        close(unit)
    end subroutine write_pins

end module m_GMSH
