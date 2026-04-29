module m_sweep_order
    use m_constants
    use m_types
    implicit none
    private ! This should be public for InitialiseGeometry
    public :: InitialiseGeometry, connectivity_and_normals ! connectivity_and_normals needs to be public for testing

contains

    subroutine InitialiseGeometry(mesh, FE, QuadSn, sweep_order)
        type(t_mesh), intent(inout)             :: mesh
        type(t_finite), intent(in)              :: FE
        type(t_sn_quadrature), intent(in)       :: QuadSn
        integer, allocatable, intent(out)       :: sweep_order(:,:)

        integer :: mm, n_angles
        real(dp) :: dir_tmp(2)

        ! 1. Build element-to-element and element-to-boundary connectivity and normals
        call connectivity_and_normals(mesh, FE)

        ! 2. Precompute the global indices for upwind neighbors to accelerate the sweep
        call precompute_upwind_indices(mesh, FE)

        ! 3. For each discrete angle, determine the element processing order (sweep order)
        n_angles = QuadSn%n_angles
        allocate(sweep_order(mesh%n_elems, n_angles))

        do mm = 1, n_angles
            dir_tmp = QuadSn%dirs(mm, 1:2)
            call generate_sweep_order(mesh, dir_tmp, sweep_order(:, mm))
        end do
    end subroutine InitialiseGeometry

    subroutine connectivity_and_normals(mesh, FE)
        type(t_mesh), intent(inout)  :: mesh
        type(t_finite), intent(in)   :: FE

        integer :: e1, e2, f1, f2, n_nodes_per_face, orientation
        logical :: found_neighbor

        n_nodes_per_face = FE%order + 1
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
                
                found_neighbor = .false.
                do e2 = e1 + 1, mesh%n_elems
                    do f2 = 1, mesh%n_faces_per_elem
                        if (mesh%face_connectivity(1, f2, e2) /= -1) cycle ! Already has a neighbor on this face

                        orientation = get_face_orientation(mesh, &
                            mesh%elems(e1, FE%face_node_map(:, f1)), &
                            mesh%elems(e2, FE%face_node_map(:, f2)))

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
            integer :: i_edge, nodes_face(FE%n_nodes_per_face)
            
            do e1 = 1, mesh%n_elems
                do f1 = 1, mesh%n_faces_per_elem
                    if (mesh%face_connectivity(1, f1, e1) == -1) then
                        nodes_face = mesh%elems(e1, FE%face_node_map(:, f1))
                        
                        do i_edge = 1, mesh%n_edges
                            ! Check if the file's Edge definition is a subset of the Element's Face
                            if (all_nodes_in_list(mesh, mesh%edges(i_edge, :), nodes_face)) then
                                mesh%face_connectivity(4, f1, e1) = mesh%boundary_ids(i_edge)
                                exit
                            end if
                        end do
                    end if
                end do
            end do
        end block

        ! Part 3: Enforce Absolute Consistency on Internal Faces
        ! This prevents cycles caused by floating point noise or "twisted" local indexing
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                e2 = mesh%face_connectivity(1, f1, e1)
                if (e2 > e1) then 
                    f2 = mesh%face_connectivity(2, f1, e1)
                    ! Ensure neighbor's normal is exactly opposite
                    mesh%face_normals(:, f2, e2) = -mesh%face_normals(:, f1, e1)
                end if
            end do
        end do
    end subroutine connectivity_and_normals
    
    subroutine precompute_upwind_indices(mesh, FE)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in)  :: FE
        integer :: ee, f, j_f, j_nf, nid, n_fac, orient

        allocate(mesh%upwind_idx(FE%n_nodes_per_face, mesh%n_faces_per_elem, mesh%n_elems))
        mesh%upwind_idx = 0

        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                nid = mesh%face_connectivity(1, f, ee)
                
                if (nid > 0) then ! Internal neighbor exists
                    n_fac  = mesh%face_connectivity(2, f, ee)
                    orient = mesh%face_connectivity(3, f, ee)
                    
                    do j_f = 1, FE%n_nodes_per_face
                        ! Handle face orientation (flipped or same)
                        if (orient == 1) then
                            j_nf = j_f
                        else
                            j_nf = FE%n_nodes_per_face - j_f + 1
                        end if
                        
                        ! Pre-calculate the absolute memory index in the angular flux array
                        ! Index = (Neighbor_ID - 1) * Nodes_Per_Elem + Local_Node_Index_on_Neighbor_Face
                        mesh%upwind_idx(j_f, f, ee) = (nid - 1) * FE%n_basis + FE%face_node_map(j_nf, n_fac)
                    end do
                end if
            end do
        end do
    end subroutine precompute_upwind_indices

    function all_nodes_in_list(mesh, subset, superset) result(is_subset)
        type(t_mesh), intent(in) :: mesh
        integer, intent(in) :: subset(:), superset(:)
        logical :: is_subset
        integer :: ii, j, k
        real(dp) :: x, y
        
        is_subset = .true.
        do ii = 1, size(subset)
            if (subset(ii) == 0) cycle ! Ignore zero-padding in edge connectivity
            x = mesh%nodes(subset(ii), 1); y = mesh%nodes(subset(ii), 2)
            j = 0
            do k = 1, size(superset)
                if (abs(mesh%nodes(superset(k),1) - x) < dp_EPSILON .and. &
                    abs(mesh%nodes(superset(k),2) - y) < dp_EPSILON) j = 1
            end do
            if (j == 0) then
                is_subset = .false.; return
            end if
        end do
    end function all_nodes_in_list

    subroutine generate_sweep_order(mesh, direction, sweep_order)
        type(t_mesh), intent(in)      :: mesh
        real(dp), intent(in)          :: direction(2)
        integer, intent(out), contiguous :: sweep_order(:)
        
        integer :: e1, e2, f1
        integer, allocatable :: queue(:), incoming(:)
        integer :: head, tail, sweep_idx, level_end

        allocate(queue(mesh%n_elems), incoming(mesh%n_elems))
        incoming = 0

        head = 1
        tail = 0
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                    if (dot_product(mesh%face_normals(:, f1, e1), direction) < -1e-12_dp) then
                    if (mesh%face_connectivity(1, f1, e1) > 0) incoming(e1) = incoming(e1) + 1
                end if
            end do
            if (incoming(e1) == 0) then
                tail = tail + 1
                queue(tail) = e1
            end if
        end do
        sweep_idx = 0
        level_end = tail
        
        do while (head <= tail)
            do while (head <= level_end)
                e1 = queue(head)
                head = head + 1
                
                sweep_idx = sweep_idx + 1
                sweep_order(sweep_idx) = e1

                do f1 = 1, mesh%n_faces_per_elem
                    if (dot_product(mesh%face_normals(:, f1, e1), direction) > 1e-12_dp) then
                        e2 = mesh%face_connectivity(1, f1, e1)
                        if (e2 > 0) then
                            incoming(e2) = incoming(e2) - 1
                            if (incoming(e2) == 0) then
                                tail = tail + 1
                                queue(tail) = e2
                            end if
                        end if
                    end if
                end do
            end do
            level_end = tail
        end do

        if (sweep_idx /= mesh%n_elems) then
             write(*,*) "Sweep Error: Elements processed:", sweep_idx, " Total:", mesh%n_elems
             stop "Sweep Error: Cycle detected in mesh or logic failure."
        end if
        
        deallocate(queue, incoming)
    end subroutine generate_sweep_order

    function get_face_orientation(mesh, nodes1, nodes2) result(orientation)
        type(t_mesh), intent(in) :: mesh
        integer, intent(in) :: nodes1(:), nodes2(:)
        integer :: orientation
        integer :: n, ii
        real(dp), allocatable :: x1(:), y1(:), x2(:), y2(:)
        
        orientation = 0
        n = size(nodes1)
        if (n /= size(nodes2)) return

        allocate(x1(n), y1(n), x2(n), y2(n))
        do ii = 1, n
            x1(ii) = mesh%nodes(nodes1(ii), 1); y1(ii) = mesh%nodes(nodes1(ii), 2)
            x2(ii) = mesh%nodes(nodes2(ii), 1); y2(ii) = mesh%nodes(nodes2(ii), 2)
        end do

        ! Compare physical coordinates to find match (plug and play for IGA knot spans)
        if (all(abs(x1 - x2) < dp_EPSILON .and. abs(y1 - y2) < dp_EPSILON)) then
            orientation = 1
        else if (all(abs(x1 - x2(n:1:-1)) < dp_EPSILON .and. abs(y1 - y2(n:1:-1)) < dp_EPSILON)) then
            orientation = -1
        end if

        deallocate(x1, y1, x2, y2)
    end function get_face_orientation
    
end module