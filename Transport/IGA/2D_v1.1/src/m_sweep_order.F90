module m_sweep_order
    use m_constants
    use m_types
    use m_basis, only: EvalNURBS2D
    implicit none
    private
    public :: InitialiseGeometry

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

        integer :: e1, e2, f1, f2, orientation
        real(dp) :: nx, ny, len, xu, yu, xv, yv, u_eval, v_eval
        real(dp) :: u1, u2, v1, v2
        real(dp), allocatable :: R(:), dR_dxi(:), dR_deta(:)
        logical :: found_neighbor

        allocate(R(FE%n_basis), dR_dxi(FE%n_basis), dR_deta(FE%n_basis))
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

                        orientation = get_face_orientation( &
                            mesh%elems(e1, pack(FE%face_node_map(:, f1, e1), FE%face_node_map(:, f1, e1) /= 0)), &
                            mesh%elems(e2, pack(FE%face_node_map(:, f2, e2), FE%face_node_map(:, f2, e2) /= 0)))

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
            integer :: i_edge
            integer, allocatable :: nodes_face(:)
            
            do e1 = 1, mesh%n_elems
                do f1 = 1, mesh%n_faces_per_elem
                    if (mesh%face_connectivity(1, f1, e1) == -1) then
                        nodes_face = mesh%elems(e1, pack(FE%face_node_map(:, f1, e1), FE%face_node_map(:, f1, e1) /= 0))
                        
                        do i_edge = 1, mesh%n_edges
                            ! Check if the file's Edge definition is a subset of the Element's Face
                            if (all_nodes_in_list(mesh%edges(i_edge, :), nodes_face)) then
                                mesh%face_connectivity(4, f1, e1) = mesh%boundary_ids(i_edge)
                                exit
                            end if
                        end do
                    end if
                end do
            end do
        end block

        ! Part 2: High-Fidelity Outward Normals Calculation for Digraph Resolution
        ! We evaluate tangents using NURBS basis functions at parametric face midpoints.
        !$omp parallel do default(none) &
        !$omp shared(mesh, FE) &
        !$omp private(e1, f1, u1, u2, v1, v2, u_eval, v_eval, R, dR_dxi, dR_deta, xu, yu, xv, yv, nx, ny, len) &
        !$omp schedule(static)
        do e1 = 1, mesh%n_elems
            u1 = mesh%span_range(1, e1); u2 = mesh%span_range(2, e1)
            v1 = mesh%span_range(3, e1); v2 = mesh%span_range(4, e1)

            do f1 = 1, mesh%n_faces_per_elem
                select case(f1)
                    case(1); u_eval = 0.5_dp * (u1 + u2); v_eval = v1
                    case(2); u_eval = u2;                v_eval = 0.5_dp * (v1 + v2)
                    case(3); u_eval = 0.5_dp * (u1 + u2); v_eval = v2
                    case(4); u_eval = u1;                v_eval = 0.5_dp * (v1 + v2)
                end select

                call EvalNURBS2D(FE, e1, mesh, u_eval, v_eval, R, dR_dxi, dR_deta)
                
                xu = dot_product(dR_dxi,  mesh%nodes(mesh%elems(e1, :), 1))
                yu = dot_product(dR_dxi,  mesh%nodes(mesh%elems(e1, :), 2))
                xv = dot_product(dR_deta, mesh%nodes(mesh%elems(e1, :), 1))
                yv = dot_product(dR_deta, mesh%nodes(mesh%elems(e1, :), 2))

                ! Outward normal logic based on parametric face orientation
                select case(f1)
                    case(1); nx =  yu; ny = -xu  ! Bottom
                    case(2); nx =  yv; ny = -xv  ! Right
                    case(3); nx = -yu; ny =  xu  ! Top
                    case(4); nx = -yv; ny =  xv  ! Left
                end select

                len = sqrt(nx**2 + ny**2)
                if (len > dp_EPSILON) then
                    mesh%face_normals(:, f1, e1) = [nx, ny] / len
                end if
            end do
        end do
        !$omp end parallel do
        deallocate(R, dR_dxi, dR_deta)
    end subroutine connectivity_and_normals
    
    subroutine precompute_upwind_indices(mesh, FE)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in)  :: FE
        integer :: ee, f, j_f, j_nf, nid, n_fac, orient, nx_ee, nx_nid

        allocate(mesh%upwind_idx(FE%n_nodes_per_face, mesh%n_faces_per_elem, mesh%n_elems))
        mesh%upwind_idx = 0

        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                if (f == 1 .or. f == 3) then; nx_ee = mesh%n_cp_xi(ee); else; nx_ee = mesh%n_cp_eta(ee); end if
                nid = mesh%face_connectivity(1, f, ee)
                
                if (nid > 0) then ! Internal neighbor exists
                    n_fac  = mesh%face_connectivity(2, f, ee)
                    orient = mesh%face_connectivity(3, f, ee)
                    if (n_fac == 1 .or. n_fac == 3) then; nx_nid = mesh%n_cp_xi(nid); else; nx_nid = mesh%n_cp_eta(nid); end if
                    
                    do j_f = 1, nx_ee
                        ! Handle face orientation (flipped or same)
                        if (orient == 1) then
                            j_nf = j_f
                        else
                            j_nf = nx_nid - j_f + 1
                        end if
                        
                        ! Pre-calculate the absolute memory index in the angular flux array
                        ! Index = (Neighbor_ID - 1) * Nodes_Per_Elem + Local_Node_Index_on_Neighbor_Face
                        mesh%upwind_idx(j_f, f, ee) = (nid - 1) * FE%n_basis + FE%face_node_map(j_nf, n_fac, nid)
                    end do
                end if
            end do
        end do
    end subroutine precompute_upwind_indices

    ! Helper to verify if boundary edges match knot span faces
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

    ! Kahn's Algorithm for Topological Sort of the Transport Digraph
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
        ! Initialize DAG: count incoming dependencies based on flow direction
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                if (dot_product(mesh%face_normals(:, f1, e1), direction) < -dp_EPSILON) then
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
        
        ! Process digraph
        do while (head <= tail)
            do while (head <= level_end)
                e1 = queue(head)
                head = head + 1
                
                sweep_idx = sweep_idx + 1
                sweep_order(sweep_idx) = e1

                do f1 = 1, mesh%n_faces_per_elem
                    if (dot_product(mesh%face_normals(:, f1, e1), direction) > dp_EPSILON) then
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
    
end module