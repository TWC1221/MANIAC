module m_sweep_order
    use m_constants
    use m_types, only: t_mesh, t_finite, t_sn_quadrature, timer_start, timer_stop
    implicit none
    private
    public :: InitialiseGeometry

    type :: t_face_info
        integer :: elem_id
        integer :: face_id
        integer, allocatable :: sorted_nodes(:)
        integer, allocatable :: original_nodes(:)
    end type t_face_info

contains

    subroutine InitialiseGeometry(mesh, FE, QuadSn, sweep_order)
        type(t_mesh), intent(inout)             :: mesh
        type(t_finite), intent(in)              :: FE
        type(t_sn_quadrature), intent(in)       :: QuadSn
        integer, allocatable, intent(out)       :: sweep_order(:,:)

        integer :: mm
        real(dp) :: dir_tmp(3)

        call timer_start('EVAL - Connectivity')
        ! 1. Build element-to-element and element-to-boundary connectivity and normals
        call connectivity_and_normals(mesh, FE)
        call timer_stop('EVAL - Connectivity')

        call timer_start('EVAL - Upwind Indices')
        ! 2. Precompute the global indices for upwind neighbors to accelerate the sweep
        call precompute_upwind_indices(mesh, FE)
        call timer_stop('EVAL - Upwind Indices')

        ! 3. For each discrete angle, determine the element processing order (sweep order)
        allocate(sweep_order(mesh%n_elems, QuadSn%n_angles))
        call timer_start('EVAL - Sweep Order')
        do mm = 1, QuadSn%n_angles
            dir_tmp = QuadSn%dirs(mm, :)
            call generate_sweep_order(mesh, dir_tmp, sweep_order(:, mm))
        end do
        call timer_stop('EVAL - Sweep Order')
        call timer_stop('Initialisation')

    end subroutine InitialiseGeometry

    subroutine connectivity_and_normals(mesh, FE)
        type(t_mesh), intent(inout)  :: mesh
        type(t_finite), intent(in)   :: FE

        integer :: e1, e2, f1, f2, k
        real(dp) :: len, v1(3), v2(3), normal(3), p1(3), p2(3), p3(3)
        real(dp) :: centroid(3), face_centroid(3), vec_out(3)
        type(t_face_info), allocatable :: all_faces(:)
        integer :: i, count, n_total_faces
        integer :: nodes_temp(FE%n_nodes_per_face)
        integer, allocatable :: edge_nodes_temp(:)

        mesh%n_faces_per_elem = 6

        allocate(mesh%face_connectivity(4, mesh%n_faces_per_elem, mesh%n_elems))
        allocate(mesh%face_normals(mesh%dim, mesh%n_faces_per_elem, mesh%n_elems))
        allocate(edge_nodes_temp(FE%n_nodes_per_face))

        ! Consolidated per-face metadata: dims = (4, n_faces_per_elem, n_elems)
        mesh%face_connectivity(1,:,:) = -1    ! face_connectivity(1, f, e) = neighbor element id (or -1)
        mesh%face_connectivity(2,:,:) = -1    ! face_connectivity(2, f, e) = neighbor face index
        mesh%face_connectivity(3,:,:) = 0     ! face_connectivity(3, f, e) = orientation (+1/-1/0)
        mesh%face_connectivity(4,:,:) = 0     ! face_connectivity(4, f, e) = boundary id (0 = internal)
        mesh%face_normals = 0.0_dp

        n_total_faces = mesh%n_elems * mesh%n_faces_per_elem
        allocate(all_faces(n_total_faces))

        count = 0
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                count = count + 1
                all_faces(count)%elem_id = e1
                all_faces(count)%face_id = f1
                
                ! Robust topological matching: use only the 4 corner nodes of the face
                ! regardless of the element order.
                allocate(all_faces(count)%original_nodes(4))
                allocate(all_faces(count)%sorted_nodes(4))
                
                all_faces(count)%original_nodes(1) = mesh%elems(e1, FE%face_node_map(1, f1))
                all_faces(count)%original_nodes(2) = mesh%elems(e1, FE%face_node_map(FE%order + 1, f1))
                all_faces(count)%original_nodes(3) = mesh%elems(e1, FE%face_node_map(FE%order*(FE%order+1) + 1, f1))
                all_faces(count)%original_nodes(4) = mesh%elems(e1, FE%face_node_map((FE%order+1)**2, f1))
                
                all_faces(count)%sorted_nodes   = all_faces(count)%original_nodes
                call sort_int_array(all_faces(count)%sorted_nodes)
            end do
        end do

        call sort_faces(all_faces, n_total_faces)

        do i = 1, n_total_faces - 1
            ! Check if current face matches next face
            if (compare_faces(all_faces(i), all_faces(i+1))) then
                ! We have a match
                e1 = all_faces(i)%elem_id
                f1 = all_faces(i)%face_id
                e2 = all_faces(i+1)%elem_id
                f2 = all_faces(i+1)%face_id

                mesh%face_connectivity(1, f1, e1) = e2
                mesh%face_connectivity(2, f1, e1) = f2
                mesh%face_connectivity(1, f2, e2) = e1
                mesh%face_connectivity(2, f2, e2) = f1
                mesh%face_connectivity(3, f1, e1) = 1 
                mesh%face_connectivity(3, f2, e2) = 1
            end if
        end do

        ! Cleanup
        do i = 1, n_total_faces
            if (allocated(all_faces(i)%sorted_nodes)) deallocate(all_faces(i)%sorted_nodes)
            if (allocated(all_faces(i)%original_nodes)) deallocate(all_faces(i)%original_nodes)
        end do
        deallocate(all_faces)
                
        ! Part 2: Assign Boundary IDs to faces
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                if (mesh%face_connectivity(1, f1, e1) == -1) then ! Boundary face
                     nodes_temp(1) = mesh%elems(e1, FE%face_node_map(1, f1))
                     nodes_temp(2) = mesh%elems(e1, FE%face_node_map(FE%order + 1, f1))
                     nodes_temp(3) = mesh%elems(e1, FE%face_node_map(FE%order*(FE%order+1) + 1, f1))
                     nodes_temp(4) = mesh%elems(e1, FE%face_node_map((FE%order+1)**2, f1))
                     call sort_int_array(nodes_temp(1:4))

                     do i = 1, mesh%n_edges
                          ! Compare against the 4 corners available in the mesh edges list
                          edge_nodes_temp(1:4) = mesh%edges(i, 1:4)
                          call sort_int_array(edge_nodes_temp(1:4))
                          if (all(nodes_temp(1:4) == edge_nodes_temp(1:4))) then
                               mesh%face_connectivity(4, f1, e1) = mesh%boundary_ids(i)
                               exit
                          end if
                     end do
                end if
            end do
        end do

        ! Part 3: Outward Normals Calculation
        !$omp parallel do default(none) &
        !$omp shared(mesh, FE) &
        !$omp private(e1, f1, v1, v2, normal, len, p1, p2, p3, k, centroid, face_centroid, vec_out) &
        !$omp schedule(static)
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                ! Calculate 3D Normal using cross product of two edges of the face
                ! Use three non-collinear nodes of the face to define two vectors
                p1 = mesh%nodes(mesh%elems(e1, FE%face_node_map(1, f1)), :)
                p2 = mesh%nodes(mesh%elems(e1, FE%face_node_map(2, f1)), :)
                p3 = mesh%nodes(mesh%elems(e1, FE%face_node_map(FE%order + 2, f1)), :)
                
                ! v1 = p2 - p1, v2 = p3 - p1
                v1 = p2 - p1
                v2 = p3 - p1

                ! Cross product v1 x v2
                normal(1) = v1(2)*v2(3) - v1(3)*v2(2)
                normal(2) = v1(3)*v2(1) - v1(1)*v2(3)
                normal(3) = v1(1)*v2(2) - v1(2)*v2(1)
                
                len = sqrt(dot_product(normal, normal))
                
                ! Ensure normal is pointing OUTWARD
                ! Compute element centroid
                centroid = 0.0_dp
                do k = 1, FE%n_basis
                    centroid = centroid + mesh%nodes(mesh%elems(e1, k), :)
                end do
                centroid = centroid / real(FE%n_basis, dp)

                ! Compute rough face centroid
                face_centroid = (p1 + p2 + p3) / 3.0_dp
                vec_out = face_centroid - centroid
                
                ! Flip if pointing inward
                if (dot_product(normal, vec_out) < 0.0_dp) normal = -normal

                if (len > dp_EPSILON) then
                    mesh%face_normals(:, f1, e1) = normal / len
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine connectivity_and_normals
    
    subroutine precompute_upwind_indices(mesh, FE)
        type(t_mesh), intent(inout) :: mesh
        type(t_finite), intent(in)  :: FE
        integer :: ee, f, j_f, j_nf, nid, n_fac
        integer :: nodes_self(FE%n_nodes_per_face), nodes_neigh(FE%n_nodes_per_face)
        integer :: glob_self, glob_neigh

        allocate(mesh%upwind_idx(FE%n_nodes_per_face, mesh%n_faces_per_elem, mesh%n_elems))
        mesh%upwind_idx = 0

        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                nid = mesh%face_connectivity(1, f, ee)
                
                if (nid > 0) then ! Internal neighbor exists
                    n_fac  = mesh%face_connectivity(2, f, ee)
                    
                    nodes_self  = mesh%elems(ee, FE%face_node_map(:, f))
                    nodes_neigh = mesh%elems(nid, FE%face_node_map(:, n_fac))
                    
                    do j_f = 1, FE%n_nodes_per_face
                        glob_self = nodes_self(j_f)
                        ! Robust search for the matching node on neighbor face
                        do j_nf = 1, FE%n_nodes_per_face
                            glob_neigh = nodes_neigh(j_nf)
                            ! Robust search: fallback to coordinate comparison if node IDs are not shared (common in 3D HO meshes)
                            if (glob_neigh == glob_self .or. &
                                all(abs(mesh%nodes(glob_self, :) - mesh%nodes(glob_neigh, :)) < 1e-8_dp)) then
                                mesh%upwind_idx(j_f, f, ee) = (nid - 1) * FE%n_basis + FE%face_node_map(j_nf, n_fac)
                                exit
                            end if
                        end do
                    end do
                end if
            end do
        end do
    end subroutine precompute_upwind_indices

    subroutine generate_sweep_order(mesh, direction, sweep_order)
        type(t_mesh), intent(in)      :: mesh
        real(dp), intent(in)          :: direction(3)
        integer, intent(out), contiguous :: sweep_order(:)
        
        integer :: e1, e2, f1
        integer, allocatable :: queue(:), incoming(:)
        integer :: head, tail, sweep_idx, level_end, i_loop

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
        
        ! Safe-guard for infinite loops
        i_loop = 0
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
            i_loop = i_loop + 1
            if (i_loop > mesh%n_elems + 100) then
                 ! Cycle detected or graph issue
                 exit 
            end if
        end do

        if (sweep_idx /= mesh%n_elems) then
             write(*,*) "Sweep Error: Elements processed:", sweep_idx, " Total:", mesh%n_elems
             stop "Sweep Error: Cycle detected in mesh or logic failure."
        end if
        
        deallocate(queue, incoming)
    end subroutine generate_sweep_order

    ! --- Sorting Utilities ---
    subroutine sort_int_array(arr)
        integer, intent(inout) :: arr(:)
        integer :: i, j, temp
        do i = 1, size(arr)-1
            do j = i+1, size(arr)
                if (arr(j) < arr(i)) then
                    temp = arr(i); arr(i) = arr(j); arr(j) = temp
                end if
            end do
        end do
    end subroutine sort_int_array

    subroutine sort_faces(faces, n)
        type(t_face_info), intent(inout) :: faces(:)
        integer, intent(in) :: n
        call quicksort_faces(faces, 1, n)
    end subroutine sort_faces

    recursive subroutine quicksort_faces(a, first, last)
        type(t_face_info), intent(inout) :: a(:)
        integer, intent(in) :: first, last
        integer :: i, j
        type(t_face_info) :: temp, pivot
        
        if (first < last) then
            pivot = a(first)
            i = first
            j = last
            do while (i < j)
                do while (compare_faces_lt(a(i), pivot) .or. compare_faces(a(i), pivot))
                    if (i >= last) exit
                    i = i + 1
                end do
                do while (compare_faces_lt(pivot, a(j)))
                    j = j - 1
                end do
                if (i < j) then
                    temp = a(i); a(i) = a(j); a(j) = temp
                end if
            end do
            temp = a(first); a(first) = a(j); a(j) = temp
            call quicksort_faces(a, first, j-1)
            call quicksort_faces(a, j+1, last)
        end if
    end subroutine quicksort_faces

    function compare_faces(f1, f2) result(res)
        type(t_face_info), intent(in) :: f1, f2
        logical :: res
        res = all(f1%sorted_nodes == f2%sorted_nodes)
    end function compare_faces

    function compare_faces_lt(f1, f2) result(res)
        type(t_face_info), intent(in) :: f1, f2
        logical :: res
        integer :: i
        res = .false.
        do i = 1, size(f1%sorted_nodes)
            if (f1%sorted_nodes(i) < f2%sorted_nodes(i)) then; res = .true.; return; end if
            if (f1%sorted_nodes(i) > f2%sorted_nodes(i)) then; res = .false.; return; end if
        end do
    end function compare_faces_lt
    
end module