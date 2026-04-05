module m_utilities
    use m_constants
    use m_types, only: t_mesh, t_finite, t_sn_quadrature, t_timer, timers, num_timers, &
                       timer_start, timer_stop, int_to_str
    use m_finite_elements, only: FE_basis_1D, GetArbitraryBasis, GetMapping
    use m_material
    use omp_lib
    implicit none
    private
    public :: print_splash, check_mesh_quality, print_timing_report, export_vtk !PinPowerNormalisation

    ! --- Mesh Quality Print Out Parameters --- 
    character(len=*), parameter :: SEP = " |-------------------------------------------------------------------------|"
    character(len=*), parameter :: HDR = " |=========================================================================|"

contains

    ! =========================================================================================================
    ! LOGO BANNER PRINT functions - just for fun and to give the code some personality
    ! =========================================================================================================

subroutine print_splash()
    implicit none
    integer :: i
    character(len=148), dimension(9) :: logo
    character(len=8)  :: date_str
    character(len=10) :: time_str

    logo(1) = "__/\\\\____________/\\\\___________/\\\\\\\\\___________/\\\\\_____/\\\________/\\\\\\\\\\\___________/\\\\\\\\\_________________/\\\\\\\\\_        "
    logo(2) = " _\/\\\\\\________/\\\\\\_________/\\\\\\\\\\\\\________\/\\\\\\___\/\\\_______\/////\\\///__________/\\\\\\\\\\\\\____________/\\\////////__       "
    logo(3) = "  _\/\\\//\\\____/\\\//\\\________/\\\/////////\\\_______\/\\\/\\\__\/\\\___________\/\\\____________/\\\/////////\\\_________/\\\/___________      "
    logo(4) = "   _\/\\\\///\\\/\\\/_\/\\\_______\/\\\_______\/\\\_______\/\\\//\\\_\/\\\___________\/\\\___________\/\\\_______\/\\\________/\\\_____________     "
    logo(5) = "    _\/\\\__\///\\\/___\/\\\_______\/\\\\\\\\\\\\\\\_______\/\\\\//\\\\/\\\___________\/\\\___________\/\\\\\\\\\\\\\\\_______\/\\\_____________    "
    logo(6) = "     _\/\\\____\///_____\/\\\_______\/\\\/////////\\\_______\/\\\_\//\\\/\\\___________\/\\\___________\/\\\/////////\\\_______\//\\\____________   "
    logo(7) = "      _\/\\\_____________\/\\\_______\/\\\_______\/\\\_______\/\\\__\//\\\\\\___________\/\\\___________\/\\\_______\/\\\________\///\\\__________  "
    logo(8) = "       _\/\\\_____________\/\\\__/\\\_\/\\\_______\/\\\__/\\\_\/\\\___\//\\\\\__/\\\__/\\\\\\\\\\\__/\\\_\/\\\_______\/\\\__/\\\____\////\\\\\\\\\_ "
    logo(9) = "        _\///______________\///__\///__\///________\///__\///__\///_____\/////__\///__\///////////__\///__\///________\///__\///________\/////////__"

    call date_and_time(date_str, time_str)

    write(*, *) ""
    do i = 1, 9
        write(*, '(A)') logo(i)
    end do

    write(*, *) ""
    write(*, '(A)') "           M U L T I G R O U P   A N G U L A R   N E U T R O N I C S : I N T E G R O - D I F F E R E N T I A L   A N A L Y T I C   C O M P U T A N T"
    write(*, '(A)') " ----------------------------------------------------------------------"
    
    ! --- SYSTEM METADATA ---
    write(*, '(A, A4, "-", A2, "-", A2, A, A2, ":", A2, ":", A2)') &
        " [ SYSTEM TIMESTAMP ] : ", date_str(1:4), date_str(5:6), date_str(7:8), &
        " @ ", time_str(1:2), time_str(3:4), time_str(5:6)
    
    write(*, '(A)') " [ OPERATOR TYPE    ] : DETERMINISTIC BOLTZMANN TRANSPORT"
    write(*, '(A)') " [ SPATIAL BASIS    ] : 3D FINITE ELEMENT METHOD (FEM)"
    write(*, '(A)') " [ ANGULAR BASIS    ] : DISCRETE ORDINATES (SN) QUADRATURE"
    write(*, '(A)') " [ ENERGY MODEL     ] : MULTIGROUP SPECTRAL DISCRETIZATION"
    write(*, '(A)') " ----------------------------------------------------------------------"
    write(*, *) ""

end subroutine print_splash

    ! =========================================================================================================
    ! Mesh Quality Analysis Functions
    ! =========================================================================================================

    subroutine check_mesh_quality(mesh, FE, printout, file_unit)
        type(t_mesh), intent(in)      :: mesh
        type(t_finite), intent(in)    :: FE
        logical, intent(in)           :: printout  ! Controls terminal output
        integer, intent(in)           :: file_unit ! Always writes to this unit
        integer  :: el, bad_jac_count, high_skew_count, u, n_out, u_list(2), p
        integer  :: degenerate_count, duplicate_node_count, n_hex, n_nodes, c(8)
        real(dp) :: min_jac, max_skew, max_aspect, el_vol
        real(dp) :: g_min_jac, g_max_skew, g_max_aspect, total_vol
        real(dp) :: avg_jac, avg_skew

        u_list(1) = file_unit
        n_out = 1
        if (printout) then
            u_list(2) = 6
            n_out = 2
        end if

        ! --- Initialization ---
        n_nodes           = mesh%n_nodes
        g_min_jac         = huge(1.0_dp)
        g_max_skew        = 0.0_dp
        g_max_aspect      = 0.0_dp
        total_vol         = 0.0_dp
        avg_jac           = 0.0_dp
        avg_skew          = 0.0_dp
        bad_jac_count     = 0
        high_skew_count   = 0
        degenerate_count  = 0

        call timer_start('CHECK - Mesh Audit')

        duplicate_node_count = check_coincident_nodes(mesh)

        n_hex = mesh%n_elems
        
        p = FE%order
        ! Identify tensor-product indices for the 8 corners of the element
        c(1) = 1
        c(2) = p + 1
        c(3) = (p+1)*p + 1
        c(4) = (p+1)**2
        c(5) = 1 + p*(p+1)**2
        c(6) = (p+1) + p*(p+1)**2
        c(7) = (p+1)*p + 1 + p*(p+1)**2
        c(8) = (p+1)**3

        do el = 1, mesh%n_elems
            if (is_degenerate(mesh%elems(el, :))) degenerate_count = degenerate_count + 1

            ! Extract corners in VTK CCW order: 0,1,2,3(bottom) 4,5,6,7(top)
            ! Mapping from Tensor to VTK corners: 1, 2, 4, 3, 5, 6, 8, 7
            min_jac = audit_hex8(mesh%elems(el, [c(1),c(2),c(4),c(3),c(5),c(6),c(8),c(7)]), &
                                 mesh, max_skew, max_aspect, el_vol)

            g_min_jac    = min(g_min_jac, min_jac)
            g_max_skew   = max(g_max_skew, max_skew)
            g_max_aspect = max(g_max_aspect, max_aspect)
            total_vol    = total_vol + el_vol
            avg_jac      = avg_jac + min_jac
            avg_skew     = avg_skew + max_skew

            if (min_jac <= 0.0_dp) bad_jac_count = bad_jac_count + 1
            if (max_skew > 0.85_dp) high_skew_count = high_skew_count + 1
        end do
        
        avg_jac = avg_jac / max(1, mesh%n_elems)
        avg_skew = avg_skew / max(1, mesh%n_elems)

        do u = 1, n_out
            write(u_list(u),*) ""

            write(u_list(u), '(A)') HDR
            write(u_list(u), '(A)') " |                       MESH QUALITY REPORT                               |"
            write(u_list(u), '(A)') HDR
            
            write(u_list(u), '(A)') " |  [1] GEOMETRY & TOPOLOGY INVENTORY                                      |"
            write(u_list(u), '(A, I15, T76, A)') " |      - Total Nodes:      ", n_nodes, "|"
            write(u_list(u), '(A, I15, T76, A)') " |      - Total Elements:   ", mesh%n_elems, "|"
            write(u_list(u), '(A, I15, T76, A)') " |        -> Hex-8 (3D):    ", n_hex, "|"
            write(u_list(u), '(A, I15, T76, A)') " |      - Duplicate Nodes:  ", duplicate_node_count, "|"
            
            
            write(u_list(u), '(A)') SEP
            write(u_list(u), '(A)') " |  [2] ELEMENT QUALITY METRICS                                            |"
            write(u_list(u), '(A, ES15.4, T76, A)') " |      - Min Jacobian:     ", g_min_jac, "|"
            write(u_list(u), '(A, ES15.4, T76, A)') " |      - Avg Jacobian:     ", avg_jac, "|"
            write(u_list(u), '(A, F15.4, T76, A)')  " |      - Max Skewness:     ", g_max_skew, "|"
            write(u_list(u), '(A, F15.4, T76, A)')  " |      - Max Aspect Ratio: ", g_max_aspect, "|"
            write(u_list(u), '(A, ES15.4, T76, A)') " |      - Computed Volume:  ", total_vol, "|"
            
            write(u_list(u), '(A)') SEP
            write(u_list(u), '(A)') " |  [3] MATERIAL ID INVENTORY                                              |"
            call report_mat_stats(mesh%material_ids, u_list(u))
            write(u_list(u), '(A)') SEP
            write(u_list(u), '(A)') " |  [4] BOUNDARY ID INVENTORY                                              |"
            call report_bc_stats(mesh%boundary_ids, u_list(u))

            write(u_list(u), '(A)') SEP
            write(u_list(u), '(A)') " |  [5] HEALTH SUMMARY                                                     |"
            if (degenerate_count == 0 .and. bad_jac_count == 0 .and. duplicate_node_count == 0) then
                write(u_list(u), '(A)') " |      >> STATUS: [ PASS ] Mesh is geometrically and topologically sound. |"
            else
                write(u_list(u), '(A)') " |      >> STATUS: [ WARNING ] Mesh issues detected!                       |"
                if (degenerate_count > 0)  write(u_list(u), '(A, I8, T76, A)') " |         - DEGENERATE:    ", degenerate_count, "|"
                if (bad_jac_count > 0)     write(u_list(u), '(A, I8, T76, A)') " |         - INVERTED (J<0):", bad_jac_count, "|"
                if (high_skew_count > 0)   write(u_list(u), '(A, I8, T76, A)') " |         - HIGH SKEW (>0.85):", high_skew_count, "|"
                if (duplicate_node_count > 0) write(u_list(u), '(A, I8, T76, A)') " |         - COINCIDENT NODES:", duplicate_node_count, "|"
            end if
            write(u_list(u), '(A)') HDR
        end do
        call timer_stop('CHECK - Mesh Audit')

    end subroutine check_mesh_quality

    function is_degenerate(ids) result(failed)
        integer, intent(in) :: ids(:)
        logical :: failed
        integer :: i, j
        failed = .false.
        do i = 1, size(ids)-1
            do j = i+1, size(ids)
                if (ids(i) == ids(j) .and. ids(i) /= 0) then
                    failed = .true.
                    return
                end if
            end do
        end do
    end function is_degenerate

    function check_coincident_nodes(mesh) result(duplicates)
        type(t_mesh), intent(in) :: mesh
        integer              :: duplicates, i, n
        real(dp), allocatable :: sums(:)
        integer, allocatable :: idx(:)
        real(dp), parameter  :: tol = 1.0d-8

        n = mesh%n_nodes
        duplicates = 0
        allocate(sums(n), idx(n))

        ! Optimization: Sorting by coordinate sum is a heuristic to avoid O(N^2).
        ! While it can have false positives (handled by the 'all' check), it identifies
        ! unshared face nodes effectively in 3D high-order meshes.
        ! NOTE: These duplicates are expected if m_GMSH is generating face nodes 
        ! per-element instead of using a global face-sharing algorithm.
        !$omp parallel do
        do i = 1, n
            sums(i) = mesh%nodes(i,1) + mesh%nodes(i,2) + mesh%nodes(i,3)
            idx(i) = i
        end do
        
        do i = 1, n
            sums(i) = mesh%nodes(i,1) + mesh%nodes(i,2) + mesh%nodes(i,3)
            idx(i) = i
        end do

        call sort_nodes(sums, idx, 1, n)

        do i = 1, n-1
            if (abs(sums(i+1) - sums(i)) < tol) then
                if (all(abs(mesh%nodes(idx(i),:) - mesh%nodes(idx(i+1),:)) < tol)) then
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

    function audit_hex8(ids, mesh, max_skew, max_aspect, vol) result(min_jac)
        integer, intent(in)  :: ids(8)
        type(t_mesh), intent(in) :: mesh
        real(dp), intent(out) :: max_skew, max_aspect, vol
        real(dp)             :: min_jac, jac_det, jacobian(3,3), dN(8,3)
        real(dp)             :: node_coords(3,8), xi, et, ze, L(3)
        integer              :: i
        
        do i = 1, 8
            node_coords(1, i) = mesh%nodes(ids(i), 1)
            node_coords(2, i) = mesh%nodes(ids(i), 2)
            node_coords(3, i) = mesh%nodes(ids(i), 3)
        end do

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

    subroutine report_bc_stats(bc_ids, unit)
        integer, intent(in) :: bc_ids(:), unit
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
            write(unit, '(A, T76, A)') " |      - No Boundary IDs (>=101) found.", "|"
        end if

        do i = 1, unique_count
            write(unit, '(A, I8, A, I10, T57, A)') " |      - BC ID       ", u_ids(i), " : ", counts(i), " boundary elements |"
        end do
        deallocate(u_ids, counts)
    end subroutine report_bc_stats

    subroutine report_mat_stats(mat_ids, file_unit)
        integer, intent(in) :: mat_ids(:), file_unit
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
            write(file_unit, '(A, I8, A, I10, T66, A)') " |      - Material ID ", u_ids(i), " : ", counts(i), " elements |"
        end do
        deallocate(u_ids, counts)
    end subroutine report_mat_stats

    ! =========================================================================================================
    ! Timing Functions - Wall Clock
    ! =========================================================================================================

    subroutine print_timing_report(printout, file_unit)
        logical, intent(in) :: printout
        integer :: i, u, u_list(2), n_out
        character(len=1024) :: display_name
        integer, intent(in) :: file_unit

        real(dp) :: total_time, avg_time, pct
        
        total_time = 0.0_dp
        
        ! Auto-stop total execution timer if it's still running
        call timer_stop('Total Execution')

        do i = 1, num_timers
            if (trim(timers(i)%name) == 'Total Execution') then
                total_time = timers(i)%accumulated_t
                exit
            end if
        end do
        if (total_time < 1.0e-9_dp) total_time = 1.0_dp
        
        u_list(1) = file_unit
        n_out = 1
        if (printout) then
            u_list(2) = 6 ! Standard output
            n_out = 2
        end if

        do u = 1, n_out

            write(u_list(u),*)
            write(u_list(u),'(A)') " |==========================================================================================|"
            write(u_list(u),'(A)') " |                                TIMING PERFORMANCE REPORT                                 |"
            write(u_list(u),'(A)') " |==========================================================================================|"
            write(u_list(u),'(A2, T5, A, T34, A8, T44, A12, T61, A12, T80, A10, T93, A1)') &
                " |", "Section", "Calls", "Time (s)", "Avg (s)", "Total (%)", "|"
            write(u_list(u),'(A)') " |==========================================================================================|"
            
            do i = 1, num_timers
                avg_time = 0.0_dp
                if (timers(i)%calls > 0) avg_time = timers(i)%accumulated_t / real(timers(i)%calls, dp)
                pct = (timers(i)%accumulated_t / total_time) * 100.0_dp
                
                select case (trim(timers(i)%name))
                case ("EVAL - Transport Sweep", "Eigen Update", "I/O  - Parse Mesh", "EVAL - Quadrature", "EVAL - Basis Functions", "I/O  - Parse Materials", &
                    "EVAL - Connectivity", "EVAL - Upwind Indices", "EVAL - Sweep Order", "EVAL - Transport Integrals", "EVAL - Reflective BC Mapping", "EVAL - LU Decomposition", &
                    "EVAL - Source")
                    
                    display_name = "  >> " // trim(timers(i)%name)
                    
                    write(u_list(u),'(A2, T5, A30, T33, I8, T43, F12.4, T60, F12.6, T78, A7, F6.2, A1, T93, A1)') &
                        " | ", display_name, timers(i)%calls, timers(i)%accumulated_t, avg_time, &
                        "   >>", pct, "%", "|"
                
                case default
                    display_name = trim(timers(i)%name)
                    
                    write(u_list(u),'(A2, T5, A30, T33, I8, T43, F12.4, T60, F12.6, T79, F9.2, A1, T93, A1)') &
                        " | ", display_name, timers(i)%calls, timers(i)%accumulated_t, avg_time, &
                        pct, "%", "|"
                end select
                
                if (trim(timers(i)%name) == 'Total Execution') then
                    write(u_list(u),'(A)') " |------------------------------------------------------------------------------------------|"
                end if
            end do
            
            write(u_list(u),'(A)') " |==========================================================================================|"
            write(u_list(u),*)
        end do

        ! The file_unit is managed by the calling program (main)
    end subroutine print_timing_report

    ! =========================================================================================================
    ! OUTPUT VTK File - via trilinear hexahedral interpolation
    ! =========================================================================================================    

    subroutine export_vtk(filename, mesh, FE, QuadSN, XPETSc, NGRP, PinPowers, is_SEM, is_adjoint, is_eigenvalue, &
        max_outer_iter, outer_iter, tol, report_unit, k_eff, flux_diff)
        character(len=*), intent(in) :: filename
        type(t_mesh),     intent(in) :: mesh
        type(t_finite),   intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: QuadSN
        real(dp),         intent(in) :: XPETSc(:,:), tol
        integer,          intent(in) :: NGRP, max_outer_iter
        logical,          intent(in) :: is_SEM, is_adjoint, is_eigenvalue
        real(dp), intent(in), optional :: PinPowers(0:), k_eff, flux_diff
        integer, intent(in) :: report_unit
        integer,          intent(in), optional :: outer_iter

        integer :: ee, g, i, j, k, unit_v
        integer :: nbasis, npts_elem, ncells_elem
        integer :: n_sub_nodes, n_sub_elems
        integer :: gid, cid, basep, p
        integer :: refine_level

        real(dp), allocatable :: xi_grid(:), eta_grid(:), zeta_grid(:)
        real(dp), allocatable :: N_eval(:), reordered_coords(:,:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)
        integer,  allocatable :: Cells(:,:)

        integer :: n000, n100, n110, n010
        integer :: n001, n101, n111, n011
        real(dp) :: p_pow
        character(len=1024) :: full_path, adj, sem
        integer :: slash_idx

        call timer_start('VTK Export')

        ! --- N-Order Logic ---
        p            = FE%order
        ! For p=2 (9-node), refine_level=3. For p=3 (16-node), refine_level=4.
        refine_level = p + 1  
        nbasis       = FE%n_basis
        npts_elem    = refine_level**3
        ncells_elem  = (refine_level-1)**3
        n_sub_nodes  = mesh%n_elems * npts_elem
        n_sub_elems  = mesh%n_elems * ncells_elem

        allocate(xi_grid(refine_level), eta_grid(refine_level), zeta_grid(refine_level))
        allocate(N_eval(nbasis), reordered_coords(nbasis, 3))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, NGRP))
        allocate(Cells(n_sub_elems, 8))

        ! Define interpolation coordinates in reference space [-1, 1]
        do i = 1, refine_level
            xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / real(p, dp)
            eta_grid(i) = xi_grid(i)
            zeta_grid(i) = xi_grid(i)
        end do

        ! 1. Interpolate geometry and solution to the sub-grid
        gid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*nbasis
            ! Get actual physical coordinates of the high-order nodes
            do i = 1, nbasis
                reordered_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
            end do
            
            ! Loop through zeta(k), eta (j) then xi (i) to match standard VTK ordering
            do k = 1, refine_level
                do j = 1, refine_level
                    do i = 1, refine_level
                        gid = gid + 1
                        call GetArbitraryBasis(FE, xi_grid(i), eta_grid(j), zeta_grid(k), N_eval)
                        
                        Xp(gid,1) = dot_product(N_eval, reordered_coords(:, 1))
                        Xp(gid,2) = dot_product(N_eval, reordered_coords(:, 2))
                        Xp(gid,3) = dot_product(N_eval, reordered_coords(:, 3))
                        
                        do g = 1, NGRP
                            Up(gid, g) = dot_product(N_eval, XPETSc(basep+1:basep+nbasis, g))
                        end do
                    end do
                end do
            end do
        end do

        ! 2. Generate sub-hex connectivity (Linear hexes for visualization)
        cid = 0
        do ee = 1, mesh%n_elems
            ! Start of this element's points in the global Xp array
            basep = (ee-1)*npts_elem 
            do k = 1, refine_level - 1
                do j = 1, refine_level - 1
                    do i = 1, refine_level - 1
                        ! Node indices in the gid sequence (1-based for now)
                        ! Current layer (k)
                        n000 = basep + (k-1)*refine_level**2 + (j-1)*refine_level + i
                        n100 = n000 + 1
                        n010 = n000 + refine_level
                        n110 = n010 + 1
                        
                        ! Next layer (k+1)
                        n001 = n000 + refine_level**2
                        n101 = n100 + refine_level**2
                        n011 = n010 + refine_level**2
                        n111 = n110 + refine_level**2
                        
                        cid = cid + 1
                        ! VTK_HEXAHEDRON ordering: 000, 100, 110, 010, 001, 101, 111, 011
                        Cells(cid, :) = [n000, n100, n110, n010, n001, n101, n111, n011]
                    end do
                end do
            end do
        end do

        ! 4. File Writing
        adj = ""
        sem = "_f"
        slash_idx = index(filename, '/', back=.true.)
        if (slash_idx == 0) slash_idx = index(filename, '\', back=.true.)
        
        if (is_adjoint) adj = "_adj"
        if (is_SEM) sem = "_s"

        if (present(outer_iter)) then
            write(*,*) ">>> Successful readout at iteration ", outer_iter
            
            if (slash_idx > 0) then
                full_path = trim(filename) // "/" // filename(slash_idx+1:) // trim(sem) // trim(adj) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSN%order)) // ".vtk"
            else
                full_path = trim(filename) // "/" // trim(filename) // trim(sem) // trim(adj) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSN%order)) // ".vtk"
            end if

        end if

        unit_v = 101
        open(unit=unit_v, file=trim(full_path), status='replace', action='write')
        write(unit_v, '(A)') "# vtk DataFile Version 3.0"
        write(unit_v, '(A)') "N-Order FE Solution"
        write(unit_v, '(A)') "ASCII"
        write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

        write(unit_v, '(A, I0, A)') "POINTS ", n_sub_nodes, " double"
        do gid = 1, n_sub_nodes
            write(unit_v, '(3F18.10)') Xp(gid,:)
        end do

        write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*9
        do cid = 1, n_sub_elems
            write(unit_v, '(9I10)') 8, Cells(cid,1)-1, Cells(cid,2)-1, Cells(cid,3)-1, Cells(cid,4)-1, &
                                       Cells(cid,5)-1, Cells(cid,6)-1, Cells(cid,7)-1, Cells(cid,8)-1
        end do

        write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
        do i = 1, n_sub_elems
            write(unit_v, '(I2)') 12 ! VTK_HEXAHEDRON
        end do

        ! Material & Pin IDs
        write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
        write(unit_v, '(A)') "SCALARS Material_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do i = 1, ncells_elem
                write(unit_v, '(I10)') mesh%material_ids(ee)
            end do
        end do
              
        write(unit_v, '(A)') "SCALARS PinID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do i = 1, ncells_elem
                write(unit_v, '(I10)') mesh%pin_ids(ee)
            end do
        end do

        if (present(PinPowers)) then
            write(unit_v, '(A)') "SCALARS Relative_Pin_Power double 1"
            write(unit_v, '(A)') "LOOKUP_TABLE default"
            do ee = 1, mesh%n_elems
                p_pow = 0.0_dp
                if (mesh%pin_ids(ee) > 0 .and. mesh%pin_ids(ee) <= ubound(PinPowers,1)) then
                    p_pow = PinPowers(mesh%pin_ids(ee))
                end if
                do i = 1, ncells_elem
                    write(unit_v, '(F18.10)') p_pow
                end do
            end do
        end if

        ! Solution Data
        write(unit_v, '(A, I10)') "POINT_DATA ", n_sub_nodes
        do g = 1, NGRP
            write(unit_v, '(A, I0)') "SCALARS Flux_Group_", g
            write(unit_v, '(A)') "double 1"
            write(unit_v, '(A)') "LOOKUP_TABLE default"
            do gid = 1, n_sub_nodes
                write(unit_v, '(F18.10)') Up(gid,g)
            end do
        end do

        close(unit_v)
        deallocate(xi_grid, eta_grid, zeta_grid, N_eval, reordered_coords, Xp, Up, Cells)
        call timer_stop('VTK Export')

        ! Write solution details to performance report if available
        if (present(outer_iter) .and. present(k_eff) .and. present(flux_diff)) then
            write(report_unit, '(A)') " "
            write(report_unit, '(A)') " |==========================================================================================|"
            write(report_unit, '(A)') " |                                 TRANSPORT SOLVE REPORT                                   |"
            write(report_unit, '(A)') " |==========================================================================================|"
            write(report_unit, '(A, T5, A, T40, A, T93, A)') " |", "Problem Type:", merge("  Eigenvalue", "Fixed Source", is_eigenvalue), "|"
            write(report_unit,'(A, T5, A, T40, ES12.4, T93, A)') " |", "Convergence Tolerance:", tol, "|"
            write(report_unit, '(A, T5, A, T35, I10, T93, A)') " |", "Iteration Threshold:", max_outer_iter, "|"
            write(report_unit, '(A, T5, A, T35, I10, T93, A)') " |", "Final Iteration:", outer_iter, "|"
            
            if (present(k_eff)) then
                write(report_unit, '(A, T5, A, T39, F15.10, T93, A)') " |", "K-Effective:", k_eff, "|"
            end if

            write(report_unit, '(A)') " |==========================================================================================|"
            write(report_unit, '(A)') " "

        end if

    end subroutine export_vtk

end module m_utilities