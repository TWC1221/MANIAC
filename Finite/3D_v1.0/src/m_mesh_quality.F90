module m_mesh_quality
    use m_constants, only: dp
    implicit none
    private
    public :: check_mesh_quality

    character(len=*), parameter :: SEP = " |-------------------------------------------------------------------------|"
    character(len=*), parameter :: HDR = " |=========================================================================|"

contains

    subroutine check_mesh_quality(nodes, elements_nodes, elements_vtk_type, num_elements, bc_ids, printout)
        real(dp), intent(in)          :: nodes(:, :)          
        integer, intent(in)           :: elements_nodes(:, :)  
        integer, intent(in)           :: elements_vtk_type(:)
        integer, intent(in)           :: num_elements
        integer, intent(in), optional :: bc_ids(:)
        logical, intent(in), optional :: printout

        integer  :: el, vtk_type, bad_jac_count, high_skew_count
        integer  :: degenerate_count, duplicate_node_count
        integer  :: n_hex, n_quad, n_other, n_nodes
        real(dp) :: min_jac, max_skew, max_aspect, el_vol
        real(dp) :: g_min_jac, g_max_skew, g_max_aspect, total_vol
        real(dp) :: avg_jac, avg_skew
        logical  :: log_active

        log_active = .false.
        if (present(printout)) log_active = printout
        if (.not. log_active) return

        ! --- Initialization ---
        n_nodes           = size(nodes, 2)
        g_min_jac         = huge(1.0_dp)
        g_max_skew        = 0.0_dp
        g_max_aspect      = 0.0_dp
        total_vol         = 0.0_dp
        avg_jac           = 0.0_dp
        avg_skew          = 0.0_dp
        bad_jac_count     = 0
        high_skew_count   = 0
        degenerate_count  = 0
        n_hex = 0; n_quad = 0; n_other = 0

        ! 1. Spatial check: Coincident nodes
        duplicate_node_count = check_coincident_nodes(nodes)

        ! 2. Element Loop: Quality and Topology
        do el = 1, num_elements
            vtk_type = elements_vtk_type(el)
            
            if (is_degenerate(elements_nodes(:, el), vtk_type)) then
                degenerate_count = degenerate_count + 1
            end if

            select case (vtk_type)
                case (12) ! VTK_HEXAHEDRON
                    n_hex = n_hex + 1
                    min_jac = audit_hex8(elements_nodes(1:8, el), nodes, max_skew, max_aspect, el_vol)
                case (9)  ! VTK_QUAD
                    n_quad = n_quad + 1
                    min_jac = audit_quad4(elements_nodes(1:4, el), nodes, max_skew, max_aspect, el_vol)
                case default
                    n_other = n_other + 1
                    cycle
            end select

            g_min_jac    = min(g_min_jac, min_jac)
            g_max_skew   = max(g_max_skew, max_skew)
            g_max_aspect = max(g_max_aspect, max_aspect)
            total_vol    = total_vol + el_vol
            avg_jac      = avg_jac + min_jac
            avg_skew     = avg_skew + max_skew

            if (min_jac <= 1.0d-12) bad_jac_count = bad_jac_count + 1
            if (max_skew > 0.85_dp) high_skew_count = high_skew_count + 1
        end do
        
        avg_jac = avg_jac / max(1, num_elements)
        avg_skew = avg_skew / max(1, num_elements)

        write(*,*) ""
        write(*, '(A)') HDR
        write(*, '(A)') " |                       MESH QUALITY REPORT                               |"
        write(*, '(A)') HDR
        
        write(*, '(A)') " |  [1] GEOMETRY & TOPOLOGY INVENTORY                                      |"
        write(*, '(A, I15, T76, A)') " |      - Total Nodes:      ", n_nodes, "|"
        write(*, '(A, I15, T76, A)') " |      - Total Elements:   ", num_elements, "|"
        write(*, '(A, I15, T76, A)') " |        -> Hex-8 (3D):    ", n_hex, "|"
        write(*, '(A, I15, T76, A)') " |        -> Quad-4 (2D):   ", n_quad, "|"
        if (n_other > 0) write(*, '(A, I15, T76, A)') " |        -> Other/Misc:    ", n_other, "|"
        write(*, '(A, I15, T76, A)') " |      - Duplicate Nodes:  ", duplicate_node_count, "|"
        
        
        write(*, '(A)') SEP
        write(*, '(A)') " |  [2] ELEMENT QUALITY METRICS                                            |"
        write(*, '(A, F15.7, T76, A)') " |      - Min Jacobian:     ", g_min_jac, "|"
        write(*, '(A, F15.7, T76, A)') " |      - Avg Jacobian:     ", avg_jac, "|"
        write(*, '(A, F15.4, T76, A)') " |      - Max Skewness:     ", g_max_skew, "|"
        write(*, '(A, F15.4, T76, A)') " |      - Max Aspect Ratio: ", g_max_aspect, "|"
        write(*, '(A, F15.4, T76, A)') " |      - Computed Volume:  ", total_vol, "|"
        
        if (present(bc_ids)) then
            write(*, '(A)') SEP
            write(*, '(A)') " |  [3] MATERIAL ID INVENTORY                                              |"
            call report_mat_stats(bc_ids)
            write(*, '(A)') SEP
            write(*, '(A)') " |  [4] BOUNDARY ID INVENTORY                                              |"
            call report_bc_stats(bc_ids)
        end if

        write(*, '(A)') SEP
        write(*, '(A)') " |  [5] HEALTH SUMMARY                                                     |"
        if (degenerate_count == 0 .and. bad_jac_count == 0 .and. duplicate_node_count == 0) then
            write(*, '(A)') " |      >> STATUS: [ PASS ] Mesh is topologically sound.                   |"
        else
            write(*, '(A)') " |      >> STATUS: [ FAIL ] Critical issues detected!                      |"
            if (degenerate_count > 0)  write(*, '(A, I8, T76, A)') " |         - DEGENERATE:    ", degenerate_count, "|"
            if (bad_jac_count > 0)     write(*, '(A, I8, T76, A)') " |         - INVERTED (J<0):", bad_jac_count, "|"
            if (high_skew_count > 0)   write(*, '(A, I8, T76, A)') " |         - HIGH SKEW (>0.85):", high_skew_count, "|"
            if (duplicate_node_count > 0) write(*, '(A, I8, T76, A)') " |         - COINCIDENT NODES:", duplicate_node_count, "|"
        end if
        write(*, '(A)') HDR
        write(*,*) ""

    end subroutine check_mesh_quality

    function is_degenerate(ids, vtk_type) result(failed)
        integer, intent(in) :: ids(:), vtk_type
        logical :: failed
        integer :: i, j, n
        failed = .false.
        n = merge(8, 4, vtk_type == 12)
        do i = 1, n-1
            do j = i+1, n
                if (ids(i) == ids(j) .and. ids(i) /= 0) then
                    failed = .true.
                    return
                end if
            end do
        end do
    end function is_degenerate

    function check_coincident_nodes(nodes) result(duplicates)
        real(dp), intent(in) :: nodes(:, :)
        integer              :: duplicates, i, n
        real(dp), allocatable :: sums(:)
        integer, allocatable :: idx(:)
        real(dp), parameter  :: tol = 1.0d-8

        n = size(nodes, 2)
        duplicates = 0
        allocate(sums(n), idx(n))
        
        do i = 1, n
            sums(i) = nodes(1,i) + nodes(2,i) + nodes(3,i)
            idx(i) = i
        end do

        call sort_nodes(sums, idx, 1, n)

        do i = 1, n-1
            if (abs(sums(i+1) - sums(i)) < tol) then
                if (all(abs(nodes(:, idx(i)) - nodes(:, idx(i+1))) < tol)) then
                    duplicates = duplicates + 1
                end if
            end if
        end do
        deallocate(sums, idx)
    end function check_coincident_nodes

    recursive subroutine sort_nodes(a, idx, first, last)
        real(dp), intent(inout) :: a(:)
        integer, intent(inout)  :: idx(:)
        integer, intent(in)     :: first, last
        integer  :: i, j
        real(dp) :: x, temp_a
        integer  :: temp_idx
        i = first; j = last
        x = a((first+last)/2)
        do
            do while (a(i) < x); i = i + 1; end do
            do while (a(j) > x); j = j - 1; end do
            if (i <= j) then
                temp_a = a(i); a(i) = a(j); a(j) = temp_a
                temp_idx = idx(i); idx(i) = idx(j); idx(j) = temp_idx
                i = i + 1; j = j - 1
            end if
            if (i > j) exit
        end do
        if (first < j) call sort_nodes(a, idx, first, j)
        if (i < last)  call sort_nodes(a, idx, i, last)
    end subroutine sort_nodes

    function audit_hex8(ids, all_coords, max_skew, max_aspect, vol) result(min_jac)
        integer, intent(in)  :: ids(8)
        real(dp), intent(in) :: all_coords(:, :)
        real(dp), intent(out) :: max_skew, max_aspect, vol
        real(dp)             :: min_jac, jac_det, jacobian(3,3), dN(8,3)
        real(dp)             :: node_coords(3,8), xi, et, ze, L(3)
        integer              :: i
        do i = 1, 8; node_coords(:, i) = all_coords(:, ids(i)); end do
        min_jac = huge(1.0_dp); max_skew = 0.0_dp
        do i = 1, 8
            xi = merge(-1.0_dp, 1.0_dp, i==1 .or. i==4 .or. i==5 .or. i==8)
            et = merge(-1.0_dp, 1.0_dp, i==1 .or. i==2 .or. i==5 .or. i==6)
            ze = merge(-1.0_dp, 1.0_dp, i<=4)
            call get_hex8_dN(xi, et, ze, dN)
            jacobian = matmul(node_coords, dN)
            jac_det = det3x3(jacobian)
            min_jac = min(min_jac, jac_det)
            max_skew = max(max_skew, calculate_skew(jacobian))
        end do
        call get_hex8_dN(0.0_dp, 0.0_dp, 0.0_dp, dN)
        vol = abs(det3x3(matmul(node_coords, dN))) * 8.0_dp
        L(1) = sqrt(sum((node_coords(:,2)-node_coords(:,1))**2))
        L(2) = sqrt(sum((node_coords(:,4)-node_coords(:,1))**2))
        L(3) = sqrt(sum((node_coords(:,5)-node_coords(:,1))**2))
        max_aspect = maxval(L) / (minval(L) + 1.0d-15)
    end function audit_hex8

    function audit_quad4(ids, all_coords, max_skew, max_aspect, vol) result(min_jac)
        integer, intent(in)  :: ids(4)
        real(dp), intent(in) :: all_coords(:, :)
        real(dp), intent(out) :: max_skew, max_aspect, vol
        real(dp)             :: min_jac, node_coords(3,4), v1(3), v2(3), norm(3), L(2)
        integer              :: i
        do i = 1, 4; node_coords(:, i) = all_coords(:, ids(i)); end do
        v1 = node_coords(:,2)-node_coords(:,1); v2 = node_coords(:,4)-node_coords(:,1)
        norm = [v1(2)*v2(3)-v1(3)*v2(2), v1(3)*v2(1)-v1(1)*v2(3), v1(1)*v2(2)-v1(2)*v2(1)]
        vol = sqrt(sum(norm**2)); min_jac = vol; max_skew = 0.0_dp
        L(1) = sqrt(sum(v1**2)); L(2) = sqrt(sum(v2**2))
        max_aspect = maxval(L) / (minval(L) + 1.0d-15)
    end function audit_quad4

    function calculate_skew(J) result(skew)
        real(dp), intent(in) :: J(3,3)
        real(dp) :: skew, mags(3)
        mags = [sqrt(sum(J(:,1)**2)), sqrt(sum(J(:,2)**2)), sqrt(sum(J(:,3)**2))] + 1.0d-15
        skew = max(abs(dot_product(J(:,1), J(:,2))/(mags(1)*mags(2))), &
                   abs(dot_product(J(:,2), J(:,3))/(mags(2)*mags(3))), &
                   abs(dot_product(J(:,3), J(:,1))/(mags(3)*mags(1))))
    end function calculate_skew

    subroutine get_hex8_dN(xi, et, ze, dN)
        real(dp), intent(in)  :: xi, et, ze
        real(dp), intent(out) :: dN(8,3)
        real(dp), parameter   :: xn(8) = [-1, 1, 1, -1, -1, 1, 1, -1], &
                                 en(8) = [-1, -1, 1, 1, -1, -1, 1, 1], &
                                 zn(8) = [-1, -1, -1, -1, 1, 1, 1, 1]
        integer :: i
        do i = 1, 8
            dN(i, 1) = 0.125_dp * xn(i) * (1.0_dp + en(i)*et) * (1.0_dp + zn(i)*ze)
            dN(i, 2) = 0.125_dp * (1.0_dp + xn(i)*xi) * en(i) * (1.0_dp + zn(i)*ze)
            dN(i, 3) = 0.125_dp * (1.0_dp + xn(i)*xi) * (1.0_dp + en(i)*et) * zn(i)
        end do
    end subroutine get_hex8_dN

    function det3x3(A) result(d)
        real(dp), intent(in) :: A(3,3)
        real(dp) :: d
        d = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
            A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
            A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
    end function det3x3

    subroutine report_bc_stats(bc_ids)
        integer, intent(in) :: bc_ids(:)
        integer :: i, j, n, unique_count
        integer, allocatable :: u_ids(:), counts(:)

        n = size(bc_ids)
        if (n == 0) return
        allocate(u_ids(n), counts(n))
        unique_count = 0

        do i = 1, n
            if (bc_ids(i) < 101) cycle
            do j = 1, unique_count
                if ((u_ids(j) == bc_ids(i))) then
                    counts(j) = counts(j) + 1
                    goto 20
                end if
            end do
            unique_count = unique_count + 1
            u_ids(unique_count) = bc_ids(i)
            counts(unique_count) = 1
        20 continue
        end do

        if (unique_count == 0) then
            write(*, '(A, T76, A)') " |      - No Boundary IDs (>=101) found.", "|"
        end if

        do i = 1, unique_count
            write(*, '(A, I8, A, I10, T69, A)') " |      - BC ID       ", u_ids(i), " : ", counts(i), " items |"
        end do
        deallocate(u_ids, counts)
    end subroutine report_bc_stats

    subroutine report_mat_stats(mat_ids)
        integer, intent(in) :: mat_ids(:)
        integer :: i, j, n, unique_count
        integer, allocatable :: u_ids(:), counts(:)

        n = size(mat_ids)
        if (n == 0) return
        allocate(u_ids(n), counts(n))
        unique_count = 0

        do i = 1, n
            if (mat_ids(i) >= 100 .or. mat_ids(i) <= 0) cycle
            do j = 1, unique_count
                if ((u_ids(j) == mat_ids(i))) then
                    counts(j) = counts(j) + 1
                    goto 30
                end if
            end do
            unique_count = unique_count + 1
            u_ids(unique_count) = mat_ids(i)
            counts(unique_count) = 1
        30 continue
        end do

        do i = 1, unique_count
            write(*, '(A, I8, A, I10, T66, A)') " |      - Material ID ", u_ids(i), " : ", counts(i), " elements |"
        end do
        deallocate(u_ids, counts)
    end subroutine report_mat_stats

end module m_mesh_quality