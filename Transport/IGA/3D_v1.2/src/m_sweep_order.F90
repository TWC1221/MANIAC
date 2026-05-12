module m_sweep_order
    use m_constants
    use m_types
    use m_quadrature, only: t_quadrature
    use m_basis, only: GetMapping3D
    implicit none
    private ! This should be public for InitialiseGeometry
    public :: InitialiseGeometry, connectivity_and_normals ! connectivity_and_normals needs to be public for testing

contains

    subroutine InitialiseGeometry(mesh, FE, QuadSn, QuadFace, QuadVol, sweep_order)
        type(t_mesh), intent(inout)             :: mesh
        type(t_finite), intent(in)              :: FE
        type(t_sn_quadrature), intent(in)       :: QuadSn
        type(t_quadrature), intent(in)          :: QuadFace, QuadVol
        integer, allocatable, intent(out)       :: sweep_order(:,:)

        integer :: mm, n_angles, ee, f
        integer :: count_int, count_ext
        real(dp) :: dir_tmp(3)

        ! 1. Build element-to-element and element-to-boundary connectivity and normals
        call connectivity_and_normals(mesh, FE, QuadFace, QuadVol)
        
        ! Diagnostic Summary
        call Print_Connectivity_Summary(mesh)

        ! 2. Precompute the global indices for upwind neighbors to accelerate the sweep
        call precompute_upwind_indices(mesh, FE)

        ! 3. For each discrete angle, determine the element processing order (sweep order)
        n_angles = QuadSn%n_angles
        allocate(sweep_order(mesh%n_elems, n_angles))

        do mm = 1, n_angles
            dir_tmp = QuadSn%dirs(mm, :)
            call generate_sweep_order(mesh, dir_tmp, sweep_order(:, mm))
        end do
    end subroutine InitialiseGeometry

    subroutine connectivity_and_normals(mesh, FE, QuadFace, QuadVol)
        type(t_mesh), intent(inout)  :: mesh
        type(t_finite), intent(in)   :: FE
        type(t_quadrature), intent(in) :: QuadFace, QuadVol

        integer :: ee, e1, e2, f, f1, f2, q, s_idx, p_id, u, v, w, nid, k_iter, n_pts, p
        real(dp) :: nodes(FE%n_basis, 3), dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
        real(dp) :: detJ, R(FE%n_basis)
        real(dp) :: u1, u2, v1, v2, w1, w2, xi_f, eta_f, zeta_f, J(3,3), dA(3), s1, s2, pos(3)
        real(dp), allocatable :: centroids(:,:,:)
        logical :: found_neighbor

        mesh%n_faces_per_elem = 6

        allocate(mesh%face_connectivity(4, 6, mesh%n_elems))
        allocate(mesh%face_normals(3, 6, mesh%n_elems))
        allocate(centroids(3, 6, mesh%n_elems))
        mesh%face_normals = 0.0_dp; centroids = 0.0_dp

        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
            u1 = mesh%elem_u_min(ee); u2 = mesh%elem_u_max(ee)
            v1 = mesh%elem_v_min(ee); v2 = mesh%elem_v_max(ee)
            w1 = mesh%elem_w_min(ee); w2 = mesh%elem_w_max(ee)

            do f = 1, 6
                dA = 0.0_dp
                n_pts = QuadFace%n_points
                do q = 1, QuadFace%n_points
                    select case(f)
                        case(1); xi_f = QuadFace%xi(q); eta_f = QuadFace%eta(q); zeta_f = -1.0_dp
                        case(2); xi_f = QuadFace%xi(q); eta_f = QuadFace%eta(q); zeta_f =  1.0_dp
                        case(3); xi_f = QuadFace%xi(q); eta_f = -1.0_dp;          zeta_f = QuadFace%eta(q)
                        case(4); xi_f = QuadFace%xi(q); eta_f =  1.0_dp;          zeta_f = QuadFace%eta(q)
                        case(5); xi_f = -1.0_dp;          eta_f = QuadFace%xi(q); zeta_f = QuadFace%eta(q)
                        case(6); xi_f =  1.0_dp;          eta_f = QuadFace%xi(q); zeta_f = QuadFace%eta(q)
                    end select

                    call GetMapping3D(FE, ee, mesh, q, QuadFace, u1, u2, v1, v2, w1, w2, nodes, &
                                     dN_dx, dN_dy, dN_dz, detJ, R, &
                                     xi_custom=xi_f, eta_custom=eta_f, zeta_custom=zeta_f, J_out=J, R_out=pos)
                                     
                    centroids(:, f, ee) = centroids(:, f, ee) + pos / real(n_pts, dp)

                    select case(f)
                        case(1,2) ! Zeta constant face: j1 x j2
                            s1 = 0.5_dp*(u2-u1); s2 = 0.5_dp*(v2-v1)
                            dA(1) = dA(1) + (J(1,2)*J(2,3) - J(1,3)*J(2,2)) * s1 * s2 * QuadFace%weights(q)
                            dA(2) = dA(2) + (J(1,3)*J(2,1) - J(1,1)*J(2,3)) * s1 * s2 * QuadFace%weights(q)
                            dA(3) = dA(3) + (J(1,1)*J(2,2) - J(1,2)*J(2,1)) * s1 * s2 * QuadFace%weights(q)
                        case(3,4) ! Eta constant face: j3 x j1
                            s1 = 0.5_dp*(w2-w1); s2 = 0.5_dp*(u2-u1)
                            dA(1) = dA(1) + (J(3,2)*J(1,3) - J(3,3)*J(1,2)) * s1 * s2 * QuadFace%weights(q)
                            dA(2) = dA(2) + (J(3,3)*J(1,1) - J(3,1)*J(1,3)) * s1 * s2 * QuadFace%weights(q)
                            dA(3) = dA(3) + (J(3,1)*J(1,2) - J(3,2)*J(1,1)) * s1 * s2 * QuadFace%weights(q)
                        case(5,6) ! Xi constant face: j2 x j3
                            s1 = 0.5_dp*(v2-v1); s2 = 0.5_dp*(w2-w1)
                            dA(1) = dA(1) + (J(2,2)*J(3,3) - J(2,3)*J(3,2)) * s1 * s2 * QuadFace%weights(q)
                            dA(2) = dA(2) + (J(2,3)*J(3,1) - J(2,1)*J(3,3)) * s1 * s2 * QuadFace%weights(q)
                            dA(3) = dA(3) + (J(2,1)*J(3,2) - J(2,2)*J(3,1)) * s1 * s2 * QuadFace%weights(q)
                    end select
                end do
                
                if (f == 1 .or. f == 3 .or. f == 5) dA = -dA ! Outward
                if (norm2(dA) > dp_EPSILON) then
                    mesh%face_normals(:, f, ee) = dA / norm2(dA)
                end if
            end do
        end do

        ! Consolidated per-face metadata: dims = (4, n_faces_per_elem, n_elems)
        mesh%face_connectivity(1,:,:) = -1    ! face_connectivity(1, f, e) = neighbor element id (or -1)
        mesh%face_connectivity(2,:,:) = -1    ! face_connectivity(2, f, e) = neighbor face index
        mesh%face_connectivity(3,:,:) = 0     ! face_connectivity(3, f, e) = orientation (+1/-1/0)
        mesh%face_connectivity(4,:,:) = BC_VACUUM ! Default to vacuum boundary

        ! Part 1: Neighbor Connectivity (Fast Tensor-Product Lookup)
        do ee = 1, mesh%n_elems
            p = mesh%elem_patch_id(ee)
            u = mesh%elem_span_indices(1, ee)
            v = mesh%elem_span_indices(2, ee)
            w = mesh%elem_span_indices(3, ee)

            ! Part 1: Bidirectional Parametric Connectivity
            ! Zeta Direction
            do k_iter = w-1, 1, -1
                nid = mesh%elem_map_to_id(p, u, v, k_iter)
                if (nid > 0) then
                    mesh%face_connectivity(1, 1, ee) = nid; mesh%face_connectivity(2, 1, ee) = 2
                    mesh%face_connectivity(3, 1, ee) = 1;   mesh%face_connectivity(4, 1, ee) = 0
                    mesh%face_connectivity(1, 2, nid) = ee; mesh%face_connectivity(2, 2, nid) = 1
                    mesh%face_connectivity(3, 2, nid) = 1;  mesh%face_connectivity(4, 2, nid) = 0
                    exit
                end if
            end do
            ! Eta Direction
            do k_iter = v-1, 1, -1
                nid = mesh%elem_map_to_id(p, u, k_iter, w)
                if (nid > 0) then
                    mesh%face_connectivity(1, 3, ee) = nid; mesh%face_connectivity(2, 3, ee) = 4
                    mesh%face_connectivity(3, 3, ee) = 1;   mesh%face_connectivity(4, 3, ee) = 0
                    mesh%face_connectivity(1, 4, nid) = ee; mesh%face_connectivity(2, 4, nid) = 3
                    mesh%face_connectivity(3, 4, nid) = 1;  mesh%face_connectivity(4, 4, nid) = 0
                    exit
                end if
            end do
            ! Xi Direction
            do k_iter = u-1, 1, -1
                nid = mesh%elem_map_to_id(p, k_iter, v, w)
                if (nid > 0) then
                    mesh%face_connectivity(1, 5, ee) = nid; mesh%face_connectivity(2, 5, ee) = 6
                    mesh%face_connectivity(3, 5, ee) = 1;   mesh%face_connectivity(4, 5, ee) = 0
                    mesh%face_connectivity(1, 6, nid) = ee; mesh%face_connectivity(2, 6, nid) = 5
                    mesh%face_connectivity(3, 6, nid) = 1;  mesh%face_connectivity(4, 6, nid) = 0
                    exit
                end if
            end do
        end do

        ! Part 2: Inter-Patch Connectivity via Topological and Geometric Matching
        ! Optimization: Only check elements that actually have "open" faces (face_connectivity == -1)
        do e1 = 1, mesh%n_elems
            if (.not. any(mesh%face_connectivity(1, :, e1) == -1)) cycle
            
            do f1 = 1, 6
                if (mesh%face_connectivity(1, f1, e1) /= -1) cycle
                
                do e2 = e1 + 1, mesh%n_elems
                    ! Optimization: Only check e2 if it also has open faces
                    if (.not. any(mesh%face_connectivity(1, :, e2) == -1)) cycle
                    
                    do f2 = 1, 6
                        if (mesh%face_connectivity(1, f2, e2) /= -1) cycle

                        found_neighbor = .false.

                        ! 1. Exact Topological Match: The "Global Knot Index" equivalent.
                        ! If the global Control Point IDs match, the elements MUST connect.
                        if (all_nodes_in_list(mesh%elems(e1, FE%face_node_map(:, f1)), &
                                              mesh%elems(e2, FE%face_node_map(:, f2)))) then
                            found_neighbor = .true.
                        ! 2. Geometric Fallback: Used for non-conforming patches where CPs aren't shared.
                        else if (norm2(centroids(:,f1,e1) - centroids(:,f2,e2)) < 1e-4_dp) then
                            if (dot_product(mesh%face_normals(:,f1,e1), mesh%face_normals(:,f2,e2)) < -0.9_dp) found_neighbor = .true.
                        end if

                        if (found_neighbor) then
                                ! Enforce bidirectional symmetry immediately
                                mesh%face_connectivity(1, f1, e1) = e2; mesh%face_connectivity(2, f1, e1) = f2
                                mesh%face_connectivity(1, f2, e2) = e1; mesh%face_connectivity(2, f2, e2) = f1
                                mesh%face_connectivity(3, f1, e1) = 1;  mesh%face_connectivity(3, f2, e2) = 1
                                mesh%face_connectivity(4, f1, e1) = 0;  mesh%face_connectivity(4, f2, e2) = 0
                                exit
                        end if
                    end do
                    if (mesh%face_connectivity(1, f1, e1) /= -1) exit
                end do
            end do
        end do

        ! --- FUNDAMENTAL TRUTH: Connectivity Symmetry Diagnostic ---
        do e1 = 1, mesh%n_elems
            do f1 = 1, 6
                e2 = mesh%face_connectivity(1, f1, e1)
                if (e2 > 0) then
                    f2 = mesh%face_connectivity(2, f1, e1)
                    if (mesh%face_connectivity(1, f2, e2) /= e1) then
                        write(*,'(A,I4,A,I4,A,I4,A,I4,A,I4)') &
                            "[CRITICAL] Asymmetric Connectivity: Elem ", e1, " face ", f1, &
                            " -> Elem ", e2, " face ", f2, " BUT that face points to Elem ", &
                            mesh%face_connectivity(1, f2, e2)
                    end if
                end if
            end do
        end do

        ! [DIAGNOSTIC] Report faces that still have no internal neighbor after geometric search
        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                if (mesh%face_connectivity(1, f1, e1) == -1) then
                    write(*,*) "[DIAG] Face Elem", e1, "Face", f1, " still has no internal neighbor after geometric search. It will be assigned a boundary condition."
                end if
            end do
        end do
        deallocate(centroids)

        do e1 = 1, mesh%n_elems
            do f1 = 1, mesh%n_faces_per_elem
                if (mesh%face_connectivity(1, f1, e1) == -1) then
                    found_neighbor = .false.
                    do s_idx = 1, size(mesh%surfaces)
                        if (all_nodes_in_list(mesh%elems(e1, FE%face_node_map(:, f1)), &
                                             mesh%surfaces(s_idx)%cp_ids)) then
                            mesh%face_connectivity(4, f1, e1) = mesh%surfaces(s_idx)%bc_id
                            found_neighbor = .true.
                            exit
                        end if
                    end do
                    if (.not. found_neighbor) mesh%face_connectivity(4, f1, e1) = BC_VACUUM
                end if
            end do
        end do

    end subroutine connectivity_and_normals
    
    subroutine Print_Connectivity_Summary(mesh)
        type(t_mesh), intent(in) :: mesh
        integer :: ee, f, n_int, n_ext
        n_int = 0; n_ext = 0
        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                if (mesh%face_connectivity(1, f, ee) > 0) then
                    n_int = n_int + 1
                else
                    n_ext = n_ext + 1
                end if
            end do
        end do
        write(*,'(A)') "--------------------------------------------------------"
        write(*,'(A, I10)') "  Total Faces:       ", mesh%n_elems * mesh%n_faces_per_elem
        write(*,'(A, I10)') "  Internal (Shared): ", n_int
        write(*,'(A, I10)') "  External (BC):     ", n_ext
        write(*,'(A)') "--------------------------------------------------------"
    end subroutine Print_Connectivity_Summary

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
                else
                    ! Boundary face: upwind index is the node on the face itself (needed for reflections)
                    do j_f = 1, FE%n_nodes_per_face
                        mesh%upwind_idx(j_f, f, ee) = (ee - 1) * FE%n_basis + FE%face_node_map(j_f, f)
                    end do
                end if
            end do
        end do
    end subroutine precompute_upwind_indices

    function all_nodes_in_list(subset, superset) result(is_subset)
        integer, intent(in) :: subset(:), superset(:)
        logical :: is_subset
        integer :: ii
        is_subset = .true.
        do ii = 1, size(subset)
            if (.not. any(superset == subset(ii))) then
                is_subset = .false.
                return
            end if
        end do
    end function all_nodes_in_list

    subroutine generate_sweep_order(mesh, direction, sweep_order)
        type(t_mesh), intent(in)      :: mesh
        real(dp), intent(in)          :: direction(3)
        integer, intent(out), contiguous :: sweep_order(:)

        integer :: e1, e2, f1, nid
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
            write(*,'(A,I0,A,I0)') "Sweep Error: Elements processed: ", sweep_idx, " / ", mesh%n_elems
            write(*,'(A,3F8.4)') "Direction: ", direction
            
            write(*,*) "--- DEPENDENCY AUTOPSY ---"
            do e1 = 1, mesh%n_elems
                if (incoming(e1) > 0) then
                    write(*,'(A,I0,A,I0,A)') "Element ", e1, " is STALLED (waiting on ", incoming(e1), " neighbors):"
                    do f1 = 1, mesh%n_faces_per_elem
                        if (dot_product(mesh%face_normals(:, f1, e1), direction) < -1e-12_dp) then
                            nid = mesh%face_connectivity(1, f1, e1)
                            if (nid > 0) then
                                if (incoming(nid) > 0) then
                                    write(*,'(A,I0,A,I0,A)') "  -> Face ", f1, " waiting on Elem ", nid, " (which is also stalled)"
                                else
                                    write(*,'(A,I0,A,I0,A)') "  -> Face ", f1, " waiting on Elem ", nid, " (which was processed? Logic Error)"
                                end if
                            else
                                write(*,'(A,I0,A)') "  -> Face ", f1, " is Boundary (Inflow)"
                            end if
                        end if
                    end do
                end if
            end do
            stop "STOP Sweep Error: Cycle detected in mesh or logic failure."
        end if
        
        deallocate(queue, incoming)
    end subroutine generate_sweep_order

    function get_face_orientation(mesh, nodes1, nodes2) result(orientation)
        type(t_mesh), intent(in) :: mesh
        integer, intent(in) :: nodes1(:), nodes2(:)
        integer :: orientation
        integer :: n, ii
        real(dp), allocatable :: x1(:), y1(:), z1(:), x2(:), y2(:), z2(:)
        
        orientation = 0
        n = size(nodes1)
        if (n /= size(nodes2)) return

        allocate(x1(n), y1(n), z1(n), x2(n), y2(n), z2(n))
        do ii = 1, n
            x1(ii) = mesh%nodes(nodes1(ii), 1); y1(ii) = mesh%nodes(nodes1(ii), 2)
            z1(ii) = mesh%nodes(nodes1(ii), 3)
            x2(ii) = mesh%nodes(nodes2(ii), 1); y2(ii) = mesh%nodes(nodes2(ii), 2)
            z2(ii) = mesh%nodes(nodes2(ii), 3)
        end do

        if (all(abs(x1 - x2) < dp_EPSILON .and. abs(y1 - y2) < dp_EPSILON .and. abs(z1 - z2) < dp_EPSILON)) then
            orientation = 1
        else if (all(abs(x1 - x2(n:1:-1)) < dp_EPSILON .and. abs(y1 - y2(n:1:-1)) < dp_EPSILON .and. abs(z1 - z2(n:1:-1)) < dp_EPSILON)) then
            orientation = -1
        end if

        deallocate(x1, y1, z1, x2, y2, z2)
    end function get_face_orientation
    
end module