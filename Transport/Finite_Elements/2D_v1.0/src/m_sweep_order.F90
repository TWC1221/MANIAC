module m_sweep_order
    use m_constants
    use m_types
    implicit none
    private
    public :: generate_sweep_order_and_colors, InitiliaseConnectivityandNormals

contains
      subroutine InitiliaseConnectivityandNormals(mesh, FE)
        type(t_mesh), intent(inout)  :: mesh
        type(t_finite), intent(in)   :: FE

        integer :: e1, e2, f1, f2, n_nodes_per_face, orientation
        real(dp) :: x1, y1, x2, y2, nx, ny, len
        logical :: found_neighbor

        n_nodes_per_face = FE%order + 1
        mesh%n_faces_per_elem = 4

    allocate(mesh%face_info(4, mesh%n_faces_per_elem, mesh%n_elems))
    allocate(mesh%face_normals(mesh%dim, mesh%n_faces_per_elem, mesh%n_elems))
    ! initialize: neighbor and neighbor_face -> -1, orientation -> 0, boundary_id -> 0
    mesh%face_info(1,:,:) = -1
    mesh%face_info(2,:,:) = -1
    mesh%face_info(3,:,:) = 0
    mesh%face_info(4,:,:) = 0
    mesh%face_normals = 0.0_dp

        ! Part 1: Neighbor Connectivity
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                if (mesh%face_info(1, f1, e1) /= -1) cycle ! Already found
                
                found_neighbor = .false.
                do e2 = e1 + 1, mesh%n_elems
                    do f2 = 1, mesh%n_faces_per_elem
                        if (mesh%face_info(1, f2, e2) /= -1) cycle ! Already has a neighbor on this face

                        orientation = get_face_orientation( &
                            mesh%elems(e1, FE%face_node_map(:, f1)), &
                            mesh%elems(e2, FE%face_node_map(:, f2)))

                        if (orientation /= 0) then
                            mesh%face_info(1, f1, e1) = e2
                            mesh%face_info(2, f1, e1) = f2
                            mesh%face_info(3, f1, e1) = orientation

                            mesh%face_info(1, f2, e2) = e1
                            mesh%face_info(2, f2, e2) = f1
                            mesh%face_info(3, f2, e2) = orientation
                            
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
                    if (mesh%face_info(1, f1, e1) == -1) then
                        nodes_face = mesh%elems(e1, FE%face_node_map(:, f1))
                        
                        do i_edge = 1, mesh%n_edges
                            if (all_nodes_in_list(nodes_face, mesh%edges(i_edge, :))) then
                                mesh%face_info(4, f1, e1) = mesh%edge_mats(i_edge)
                                exit
                            end if
                        end do
                    end if
                end do
            end do
        end block

        ! Part 2: Outward Normals Calculation
        !$omp parallel do default(none) &
        !$omp shared(mesh, FE, n_nodes_per_face) &
        !$omp private(e1, f1, x1, y1, x2, y2, nx, ny, len) &
        !$omp schedule(static)
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                x1 = mesh%nodes(mesh%elems(e1, FE%face_node_map(1, f1)), 1)
                y1 = mesh%nodes(mesh%elems(e1, FE%face_node_map(1, f1)), 2)
                x2 = mesh%nodes(mesh%elems(e1, FE%face_node_map(n_nodes_per_face, f1)), 1)
                y2 = mesh%nodes(mesh%elems(e1, FE%face_node_map(n_nodes_per_face, f1)), 2)

                nx = y2 - y1
                ny = x1 - x2
                len = sqrt(nx**2 + ny**2)
                
                if (len > dp_EPSILON) then
                    mesh%face_normals(:, f1, e1) = [nx, ny] / len
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine InitiliaseConnectivityandNormals

    pure function all_nodes_in_list(subset, superset) result(is_subset)
        integer, intent(in) :: subset(:), superset(:)
        logical :: is_subset
        integer :: i
        
        is_subset = .true.
        do i = 1, size(subset)
            if (.not. any(superset == subset(i))) then
                is_subset = .false.; return
            end if
        end do
    end function

    subroutine generate_sweep_order_and_colors(mesh, direction, sweep_order, colors, num_colors)
        type(t_mesh), intent(in)      :: mesh
        real(dp), intent(in)          :: direction(2)
        integer, intent(out), contiguous :: sweep_order(:)
        integer, intent(out), contiguous :: colors(:)
        integer, intent(out)          :: num_colors
        
        integer :: e1, e2, f1
        integer, allocatable :: queue(:), incoming(:)
        integer :: head, tail, sweep_idx, level_end, color_val

        allocate(queue(mesh%n_elems), incoming(mesh%n_elems))
        incoming = 0
        colors = 0

        head = 1
        tail = 0
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                    if (dot_product(mesh%face_normals(:, f1, e1), direction) < -1e-12_dp) then
                    if (mesh%face_info(1, f1, e1) > 0) incoming(e1) = incoming(e1) + 1
                end if
            end do
            if (incoming(e1) == 0) then
                tail = tail + 1
                queue(tail) = e1
            end if
        end do

        sweep_idx = 0
        color_val = 0
        level_end = tail
        
        do while (head <= tail)
            color_val = color_val + 1
            do while (head <= level_end)
                e1 = queue(head)
                head = head + 1
                
                sweep_idx = sweep_idx + 1
                sweep_order(sweep_idx) = e1
                colors(e1) = color_val

                do f1 = 1, mesh%n_faces_per_elem
                    if (dot_product(mesh%face_normals(:, f1, e1), direction) > 1e-12_dp) then
                        e2 = mesh%face_info(1, f1, e1)
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
        num_colors = color_val

        if (sweep_idx /= mesh%n_elems) then
             write(*,*) "Sweep Error: Elements processed:", sweep_idx, " Total:", mesh%n_elems
             stop "Sweep Error: Cycle detected in mesh or logic failure."
        end if
        
        deallocate(queue, incoming)
    end subroutine generate_sweep_order_and_colors

    pure function get_face_orientation(nodes1, nodes2) result(orientation)
        integer, intent(in) :: nodes1(:), nodes2(:)
        integer :: orientation
        integer :: n_nodes, i, count
        
        orientation = 0
        n_nodes = size(nodes1)
        if (n_nodes /= size(nodes2)) return

        count = 0
        do i = 1, n_nodes
            if (any(nodes2 == nodes1(i))) count = count + 1
        end do
        if (count /= n_nodes) return

        if (nodes1(1) == nodes2(1) .and. nodes1(n_nodes) == nodes2(n_nodes)) then
            orientation = 1
            return
        end if

        if (nodes1(1) == nodes2(n_nodes) .and. nodes1(n_nodes) == nodes2(1)) then
            orientation = -1
            return
        end if
    end function get_face_orientation
    
end module