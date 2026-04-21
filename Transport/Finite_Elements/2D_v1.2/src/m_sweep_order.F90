module m_sweep_order
    use m_constants
    use m_types
    implicit none
    private
    public :: InitialiseGeometry, connectivity_and_normals, precompute_upwind_indices, generate_sweep_order, get_face_orientation, get_face_cp_indices

contains

    subroutine InitialiseGeometry(mesh, FE, QuadSn, sweep_order)
        type(t_mesh), intent(inout)             :: mesh
        type(t_finite), intent(in)              :: FE
        type(t_sn_quadrature), intent(in)       :: QuadSn
        integer, allocatable, intent(out)       :: sweep_order(:,:)

        integer :: mm, n_angles, i
        real(dp) :: dir_tmp(2)

        ! 1. Build element-to-element and element-to-boundary connectivity and normals
        call connectivity_and_normals(mesh, FE)
        write(*,*) ">>> [DEBUG] Geometry Connectivity built. Checking for orphans..."

        ! 2. Precompute the global indices for upwind neighbors to accelerate the sweep
        call precompute_upwind_indices(mesh, FE)

        ! 3. For each discrete angle, determine the element processing order (sweep order)
        !    and group elements into "colors" for parallel execution.
        n_angles = QuadSn%n_angles
        allocate(sweep_order(mesh%n_elems, n_angles))

        write(*,*) ">>> [DEBUG] Generating sweep orders for ", n_angles, " angles."
        do mm = 1, n_angles
            dir_tmp = QuadSn%dirs(mm, 1:2)
            call generate_sweep_order(mesh, dir_tmp, sweep_order(:, mm))
            
            ! Robust check: Ensure uniqueness and completeness
            block
                integer :: counts(mesh%n_elems), missing_count
                counts = 0
                missing_count = 0
                do i = 1, mesh%n_elems
                    if (sweep_order(i, mm) > 0 .and. sweep_order(i, mm) <= mesh%n_elems) then
                        counts(sweep_order(i, mm)) = counts(sweep_order(i, mm)) + 1
                    else
                        write(*,*) ">>> [DEBUG] FATAL: Angle ", mm, " pos ", i, " contains invalid element ID: ", sweep_order(i, mm)
                    end if
                end do

                do i = 1, mesh%n_elems
                    if (counts(i) == 0) then
                        if (missing_count < 10) write(*,*) ">>> [DEBUG] Element ", i, " is MISSING from sweep order for angle ", mm
                        missing_count = missing_count + 1
                    else if (counts(i) > 1) then
                        write(*,*) ">>> [DEBUG] Element ", i, " is DUPLICATED (", counts(i), " times) in sweep order for angle ", mm
                    end if
                end do

                if (missing_count > 0) then
                    write(*,*) ">>> [DEBUG] Total missing elements for angle ", mm, ": ", missing_count
                end if

                if (any(counts /= 1)) then
                    write(*,*) "FATAL: Sweep order corruption for angle ", mm, ". Missing or duplicate elements."
                    stop
                end if
            end block
        end do
    end subroutine InitialiseGeometry

    subroutine connectivity_and_normals(mesh, FE)
        type(t_mesh), intent(inout)  :: mesh
        type(t_finite), intent(in)   :: FE

        integer :: e1, e2, f1, f2, orientation, n1, n2
        real(dp) :: x1, y1, x2, y2, nx, ny, len, xm, ym
        logical :: found_neighbor
        integer :: idx1(mesh%nloc), idx2(mesh%nloc), n_f1, n_f2
        integer, allocatable :: nodes1(:), nodes2(:)

        mesh%n_faces_per_elem = 4

        allocate(mesh%face_connectivity(4, mesh%n_faces_per_elem, mesh%n_elems))
        allocate(mesh%face_normals(mesh%dim, mesh%n_faces_per_elem, mesh%n_elems))

        ! Consolidated per-face metadata: dims = (4, n_faces_per_elem, n_elems)
        mesh%face_connectivity(1,:,:) = -1    ! face_connectivity(1, f, e) = neighbor element id (or -1)
        mesh%face_connectivity(2,:,:) = -1    ! face_connectivity(2, f, e) = neighbor face index
        mesh%face_connectivity(3,:,:) = 0     ! face_connectivity(3, f, e) = orientation (+1/-1/0)
        mesh%face_connectivity(4,:,:) = 0     ! face_connectivity(4, f, e) = boundary id (0 = internal)
        mesh%face_normals = 0.0_dp

        ! Part 1: Neighbor Connectivity
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                if (mesh%face_connectivity(1, f1, e1) /= -1) cycle ! Already found
                
                call get_face_cp_indices(e1, f1, mesh, idx1, n_f1)
                if (allocated(nodes1)) deallocate(nodes1)
                allocate(nodes1(n_f1))
                nodes1 = mesh%elems(e1, idx1(1:n_f1))
                
                found_neighbor = .false.
                do e2 = e1 + 1, mesh%n_elems
                    do f2 = 1, mesh%n_faces_per_elem
                        if (mesh%face_connectivity(1, f2, e2) /= -1) cycle ! Already has a neighbor on this face

                        call get_face_cp_indices(e2, f2, mesh, idx2, n_f2)
                        if (allocated(nodes2)) deallocate(nodes2)
                        allocate(nodes2(n_f2))
                        nodes2 = mesh%elems(e2, idx2(1:n_f2))

                        orientation = get_face_orientation(nodes1, nodes2)

                        if (orientation /= 0) then
                            mesh%face_connectivity(1, f1, e1) = e2
                            mesh%face_connectivity(2, f1, e1) = f2
                            mesh%face_connectivity(3, f1, e1) = orientation

                            mesh%face_connectivity(1, f2, e2) = e1
                            mesh%face_connectivity(2, f2, e2) = f1
                            mesh%face_connectivity(3, f2, e2) = orientation
                            
                            found_neighbor = .true.
                            exit ! exit f2 loop
                        end if
                    end do
                    if (found_neighbor) exit ! exit e2 loop
                end do
            end do
        end do

        block
            integer :: i_edge, n_f
            integer, allocatable :: nodes_face(:)
            
            do e1 = 1, mesh%n_elems
                do f1 = 1, mesh%n_faces_per_elem
                    if (mesh%face_connectivity(1, f1, e1) == -1) then
                        call get_face_cp_indices(e1, f1, mesh, idx1, n_f)
                        if (allocated(nodes_face)) deallocate(nodes_face)
                        allocate(nodes_face(n_f))
                        nodes_face = mesh%elems(e1, idx1(1:n_f))
                        
                        do i_edge = 1, mesh%n_edges
                            ! Check if the file's Edge definition is a subset of the Element's Face
                            if (all_nodes_in_list(mesh%edges(i_edge, :), nodes_face)) then
                                mesh%face_connectivity(4, f1, e1) = mesh%edge_mats(i_edge)
                                exit
                            end if
                        end do
                    end if
                end do
            end do
        end block

        ! Part 2: Outward Normals Calculation
        do e1 = 1, mesh%n_elems
            ! Calculate approximate patch centroid
            x1 = 0.0_dp; y1 = 0.0_dp; n1 = 0
            do n2 = 1, mesh%nloc
                if (mesh%elems(e1, n2) > 0) then
                    x1 = x1 + mesh%nodes(mesh%elems(e1, n2), 1)
                    y1 = y1 + mesh%nodes(mesh%elems(e1, n2), 2)
                    n1 = n1 + 1
                end if
            end do
            if (n1 > 0) then; x1 = x1 / real(n1, dp); y1 = y1 / real(n1, dp); end if

            do f1 = 1, mesh%n_faces_per_elem
                ! If neighbor already calculated normal, use the negative to ensure acyclic dependencies
                e2 = mesh%face_connectivity(1, f1, e1)
                if (e2 > 0 .and. e2 < e1) then
                    f2 = mesh%face_connectivity(2, f1, e1)
                    mesh%face_normals(:, f1, e1) = -mesh%face_normals(:, f2, e2)
                    cycle
                end if

                call get_face_cp_indices(e1, f1, mesh, idx1, n_f1)
                n1 = mesh%elems(e1, idx1(1))
                n2 = mesh%elems(e1, idx1(n_f1))

                if (n1 <= 0 .or. n2 <= 0) cycle

                ! Chord vector and initial normal candidate
                x2 = mesh%nodes(n2, 1) - mesh%nodes(n1, 1)
                y2 = mesh%nodes(n2, 2) - mesh%nodes(n1, 2)
                nx = y2; ny = -x2

                ! Midpoint of face chord
                xm = 0.5_dp * (mesh%nodes(n1, 1) + mesh%nodes(n2, 1))
                ym = 0.5_dp * (mesh%nodes(n1, 2) + mesh%nodes(n2, 2))

                ! Ensure normal points outward from centroid
                if (nx * (xm - x1) + ny * (ym - y1) < 0.0_dp) then
                    nx = -nx; ny = -ny
                end if

                len = sqrt(nx**2 + ny**2)
                if (len < dp_EPSILON) len = 1.0_dp
                mesh%face_normals(:, f1, e1) = [nx, ny] / len

            end do
        end do
    end subroutine connectivity_and_normals
    
    subroutine precompute_upwind_indices(mesh, FE)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in)  :: FE
        integer :: ee, f, j_f, j_nf, nid, n_fac, orient, n_f_cp, n_nf_cp
        integer :: idx_f(mesh%nloc), idx_nf(mesh%nloc)

        allocate(mesh%upwind_idx(mesh%nloc, mesh%n_faces_per_elem, mesh%n_elems))
        mesh%upwind_idx = 0

        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                nid = mesh%face_connectivity(1, f, ee)
                
                if (nid > 0) then ! Internal neighbor exists
                    n_fac  = mesh%face_connectivity(2, f, ee)
                    orient = mesh%face_connectivity(3, f, ee)
                    
                    call get_face_cp_indices(ee, f, mesh, idx_f, n_f_cp)
                    call get_face_cp_indices(nid, n_fac, mesh, idx_nf, n_nf_cp)

                    do j_f = 1, n_f_cp
                        ! Handle face orientation (flipped or same)
                        if (orient == 1) then
                            j_nf = j_f
                        else
                            j_nf = n_nf_cp - j_f + 1
                        end if
                        
                        ! Pre-calculate the absolute memory index in the angular flux array
                        ! Index = (Neighbor_ID - 1) * Nodes_Per_Elem + Local_Node_Index_on_Neighbor_Face
                        mesh%upwind_idx(j_f, f, ee) = (nid - 1) * FE%n_basis + idx_nf(j_nf)
                    end do
                end if
            end do
        end do
    end subroutine precompute_upwind_indices

    pure function all_nodes_in_list(subset, superset) result(is_subset)
        integer, intent(in) :: subset(:), superset(:)
        logical :: is_subset
        integer :: i
        
        is_subset = .true.
        do i = 1, size(subset)
            if (subset(i) == 0) cycle ! Ignore zero-padding in edge connectivity
            if (.not. any(superset == subset(i))) then
                is_subset = .false.; return
            end if
        end do
    end function all_nodes_in_list

    subroutine generate_sweep_order(mesh, direction, sweep_order)
        type(t_mesh), intent(in)      :: mesh
        real(dp), intent(in)          :: direction(2)
        integer, intent(out)          :: sweep_order(:)
        
        integer :: i, f, neighbor, n, head, tail
        integer, allocatable :: in_degree(:), queue(:)
        real(dp) :: o_n

        n = mesh%n_elems
        allocate(in_degree(n), queue(n))
        in_degree = 0
        sweep_order = 0

        ! 1. Calculate in-degrees based on upwind neighbors for this direction
        do i = 1, n
            do f = 1, mesh%n_faces_per_elem
                neighbor = mesh%face_connectivity(1, f, i)
                if (neighbor > 0) then
                    o_n = dot_product(direction, mesh%face_normals(:, f, i))
                    ! If o_n < 0, this is an INFLOW face for element i
                    if (o_n < -1e-12_dp) in_degree(i) = in_degree(i) + 1
                end if
            end do
        end do

        ! 2. Initialize queue with elements that have no upwind neighbors
        head = 1; tail = 0
        do i = 1, n
            if (in_degree(i) == 0) then
                tail = tail + 1
                queue(tail) = i
            end if
        end do

        ! 3. Process the queue (Topological Sort)
        i = 0
        do while (head <= tail)
            i = i + 1
            sweep_order(i) = queue(head)
            
            ! For the current element, find its downwind neighbors and reduce their in-degree
            do f = 1, mesh%n_faces_per_elem
                neighbor = mesh%face_connectivity(1, f, queue(head))
                if (neighbor > 0) then
                    o_n = dot_product(direction, mesh%face_normals(:, f, queue(head)))
                    ! If o_n > 0, this is an OUTFLOW face for queue(head), hence INFLOW for neighbor
                    if (o_n > 1e-12_dp) then
                        in_degree(neighbor) = in_degree(neighbor) - 1
                        if (in_degree(neighbor) == 0) then
                            tail = tail + 1
                            queue(tail) = neighbor
                        end if
                    end if
                end if
            end do
            head = head + 1
        end do

        if (i < n) then
            write(*,*) "WARNING: Cyclic dependency detected in mesh for direction ", direction
            ! Fallback to centroid projection if topological sort stalls
            call fallback_centroid_sort(mesh, direction, sweep_order)
        end if
        
        deallocate(in_degree, queue)
    contains
        subroutine fallback_centroid_sort(m, dir, order)
            type(t_mesh), intent(in) :: m
            real(dp), intent(in) :: dir(2)
            integer, intent(out) :: order(:)
            real(dp) :: c(2), p(m%n_elems)
            integer :: ii, jj, kk, node
            do ii = 1, m%n_elems
                c = 0.0_dp; kk = 0
                do jj = 1, m%nloc
                    node = m%elems(ii, jj)
                    if (node > 0) then; c = c + m%nodes(node, 1:2); kk = kk + 1; end if
                end do
                if (kk > 0) c = c / real(kk, dp)
                p(ii) = dot_product(c, dir)
                order(ii) = ii
            end do
            ! Simple selection sort for fallback
            do ii = 1, m%n_elems - 1
                kk = ii
                do jj = ii + 1, m%n_elems
                    if (p(order(jj)) < p(order(kk))) kk = jj
                end do
                node = order(ii); order(ii) = order(kk); order(kk) = node
            end do
        end subroutine fallback_centroid_sort
    end subroutine generate_sweep_order

    pure function get_face_orientation(nodes1, nodes2) result(orientation)
        integer, intent(in) :: nodes1(:), nodes2(:)
        integer :: orientation
        integer :: n_nodes
        
        orientation = 0
        n_nodes = size(nodes1)
        if (n_nodes /= size(nodes2)) return

        ! Check for same orientation by comparing all nodes in the sequence.
        if (all(nodes1 == nodes2)) then
            orientation = 1
            return
        end if

        ! Check for flipped orientation by comparing all nodes in the reversed sequence.
        if (all(nodes1 == nodes2(n_nodes:1:-1))) then
            orientation = -1
            return
        end if

    end function get_face_orientation
    
    subroutine get_face_cp_indices(ee, f, mesh, idxs, n_f_cp)
        integer, intent(in) :: ee, f
        type(t_mesh), intent(in) :: mesh
        integer, intent(out) :: idxs(:)
        integer, intent(out) :: n_f_cp
        integer :: i, nxi, neta
        nxi = mesh%n_cp_xi(ee); neta = mesh%n_cp_eta(ee)
        select case(f)
        case(1); n_f_cp = neta; do i=1,neta; idxs(i)=(i-1)*nxi+1; end do ! Xi_min
        case(2); n_f_cp = neta; do i=1,neta; idxs(i)=i*nxi; end do       ! Xi_max
        case(3); n_f_cp = nxi;  do i=1,nxi;  idxs(i)=i; end do           ! Eta_min
        case(4); n_f_cp = nxi;  do i=1,nxi;  idxs(i)=(neta-1)*nxi+i; end do ! Eta_max
        end select
    end subroutine

end module