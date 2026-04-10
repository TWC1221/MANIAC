module m_outpVTK_dfem
    use m_constants
    use m_types
    use m_finite_elements
    use m_material
    use m_quadrature
    implicit none
    private
    public :: export_dfem_vtk, PinPowerNormalisation, derive_case_nametag, int_to_str


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
        dot_pos = index(tag, '.asmg', back=.true.)
        if (dot_pos > 0) then
            tag = tag(:dot_pos-1)
        end if
    end function derive_case_nametag
    
    subroutine PinPowerNormalisation(mesh, FE, materials, Quad2D, n_groups, &
                                     scalar_flux, angular_flux, total_src, pin_powers, filename)
        type(t_mesh), intent(inOUT)                    :: mesh
        type(t_finite), intent(in)                  :: FE
        type(t_material), intent(in)                :: materials(:)
        type(t_quadrature), intent(in)              :: Quad2D
        integer, intent(in)                         :: n_groups
        real(dp), intent(inout)                     :: scalar_flux(:,:), angular_flux(:,:,:), total_src(:,:)
        real(dp), allocatable, intent(inout)        :: pin_powers(:)
        character(len=*), intent(in)                :: filename
        real(dp), allocatable                       :: pin_powers_ref(:)

        integer :: ee, g, k, i, pin_cnt, idx_start, idx_end
        real(dp) :: fission_rate_total, elem_fission, factor, detJ
        real(dp) :: el_nodes(FE%n_basis, 2), el_weights(FE%n_basis)
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis)
        real(dp) :: u_local(n_groups), basis_rational(FE%n_basis)
        real(dp) :: total_pin_power, min_pin_power, max_err
        integer :: max_pin

        max_pin = maxval(mesh%pin_ids)
        allocate(pin_powers(0:max_pin), pin_powers_ref(0:max_pin)); pin_powers = 0.0_dp

        fission_rate_total = 0.0_dp

        do ee = 1, mesh%n_elems
            if (any(mesh%material_ids(ee) == [1,6,7,8])) mesh%pin_ids(ee) = 0 ! Assign to ID 0 to avoid overlapping with fuel pins 1...N
            if (.not. allocated(materials(mesh%material_ids(ee))%SigF)) cycle
            if (sum(materials(mesh%material_ids(ee))%SigF) < 1.0e-20) cycle

            el_nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), 1:2)
            el_weights = mesh%weights(mesh%elems(ee, 1:FE%n_basis))

            idx_start = (ee-1)*FE%n_basis + 1
            idx_end   = ee*FE%n_basis
            
            elem_fission = 0.0_dp
            do k = 1, Quad2D%n_points
                call GetMapping(FE, k, el_nodes, dN_dx, dN_dy, detJ, el_weights, basis_rational)
                if (abs(detJ) < 1.0d-15) cycle

                u_local = 0.0_dp
                do g = 1, n_groups
                    u_local(g) = dot_product(scalar_flux(idx_start:idx_end, g), basis_rational)
                end do

                elem_fission = elem_fission + Quad2D%weights(k) * abs(detJ) * dot_product(materials(mesh%material_ids(ee))%SigF, u_local)
            end do
            fission_rate_total = fission_rate_total + elem_fission
            if (mesh%pin_ids(ee) > 0) pin_powers(mesh%pin_ids(ee)) = pin_powers(mesh%pin_ids(ee)) + elem_fission
        end do

        pin_cnt = 0
        total_pin_power = 0.0_dp
        do i = 1, max_pin
            if (pin_powers(i) > 1.0e-12_dp) then
                 pin_cnt = pin_cnt + 1
                 total_pin_power = total_pin_power + pin_powers(i)
            end if
        end do

        write(*,*) ">>> Normalising to 1 fission/sec/pin for ", trim(int_to_str(pin_cnt)), " active pins."

        call ParsePinPowers(pin_powers_ref, filename)

        min_pin_power = 1.0_dp
        max_err = 0.0_dp
        if (total_pin_power > 1.0e-12_dp) then
            factor =  real(pin_cnt, dp)/total_pin_power
            scalar_flux  = scalar_flux * factor
            angular_flux = angular_flux * factor
            total_src    = total_src * factor
            pin_powers   = pin_powers * factor
            
            do i = 1, min(max_pin, ubound(pin_powers_ref, 1))
                if (pin_powers(i) >= 0.01) min_pin_power = min(min_pin_power, pin_powers(i))
                if (pin_powers_ref(i) >= 0.01) then
                    max_err = max(max_err, abs(pin_powers(i) - pin_powers_ref(i)) / pin_powers_ref(i))
                end if
            end do

            write(*,*) ">>> Max pin power:", maxval(pin_powers)
            write(*,*) ">>> Min pin power:", min_pin_power
            write(*, '(">>> Max per cent error: ", F8.4, " %")') max_err*100.0_dp
        end if
    end subroutine PinPowerNormalisation


    subroutine ParsePinPowers(reference_powers, filename)
        real(dp), allocatable, intent(out) :: reference_powers(:)
        character(len=*), intent(in) :: filename

        real(dp), allocatable :: raw(:,:)
        integer :: i, j, unit, ios, n = 34, m = 17
        integer :: Ax, Ay, i_loc, j_loc, assm_idx
        integer :: max_pin_id, pin_id

        allocate(raw(n, n))
        raw = 0.0_dp
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) "Error: Cannot open pin power table: ", trim(filename)
            if (allocated(raw)) deallocate(raw)
            return
        end if

        ! Read entire 34x34 data block. Fortran '*' read is robust to line breaks.
        read(unit, *, iostat=ios) raw
        close(unit)
            do j = 1, n
                do i = 1, n
                    raw(i, j) = max(raw(i, j), raw(j, i))
                    print*, "Read raw power for (i=", i, ", j=", j, "): ", raw(i, j) ! Debug print
                end do
            end do

        ! Map 2D grid to 1D pin_id using block-sequential ordering:
        ! Assemblies are indexed BL=1, BR=2, TL=3, TR=4
        max_pin_id = 2312
        allocate(reference_powers(0:max_pin_id))
        reference_powers = 0.0_dp

        do j = 1, n
            Ay    = (j - 1) / m + 1   ! Map j=1..17 (Top) to Ay=1; j=18..34 (Bottom) to Ay=2
            j_loc = m - mod(j - 1, m) ! Intra-assembly flip: file row 1 (Top) -> mesh j_loc 17 (Top)
            
            do i = 1, n
                Ax    = (i - 1) / m + 1
                i_loc = mod(i - 1, m) + 1 
                
                ! Map 2x2 grid to mesh assemblies 2,3,4,5 (Sequential block ordering)
                ! Using Ay=1 for top half results in assm_idx 2,3. Ay=2 for bottom half results in 4,5.
                assm_idx = (Ay - 1) * 2 + Ax + 1
                pin_id   = (assm_idx - 1) * (m * m) + (j_loc - 1) * m + i_loc + 867
                
                if (pin_id <= max_pin_id) reference_powers(pin_id) = raw(j, i)

                PRINT*  , "Pin ID:", pin_id, " (i=", i, ", j=", j, ") -> Assm:", assm_idx, " (i_loc=", i_loc, ", j_loc=", j_loc, ") with power ", raw(j,i)
            end do
        end do
        deallocate(raw)
    end subroutine ParsePinPowers

    subroutine export_dfem_vtk(filename, mesh, FE, QuadSN, XPETSc, NGRP, PinPowers, is_SEM, is_adjoint)
        character(len=*), intent(in) :: filename
        type(t_mesh),     intent(in) :: mesh
        type(t_finite),   intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: QuadSN
        real(dp),         intent(in) :: XPETSc(:,:) 
        integer,          intent(in) :: NGRP
        logical,          intent(in) :: is_SEM, is_adjoint
        real(dp), intent(in), optional :: PinPowers(0:)

        integer :: ee, g, i, j, unit_v
        integer :: nbasis, npts_elem, ncells_elem
        integer :: n_sub_nodes, n_sub_elems
        integer :: gid, cid, basep, p
        integer :: refine_level

        real(dp), allocatable :: xi_grid(:), eta_grid(:), weights_local(:)
        real(dp), allocatable :: N_eval(:), reordered_coords(:,:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)
        real(dp)              :: p_pow
        integer,  allocatable :: Cells(:,:)

        integer :: n00, n10, n11, n01
        character(len=1024) :: full_path, adj, sem
        integer :: slash_idx

        ! --- N-Order Logic ---
        p            = FE%order
        ! For p=2 (9-node), refine_level=3. For p=3 (16-node), refine_level=4.
        refine_level = p + 1  
        nbasis       = FE%n_basis
        npts_elem    = refine_level**2
        ncells_elem  = (refine_level-1)**2
        n_sub_nodes  = mesh%n_elems * npts_elem
        n_sub_elems  = mesh%n_elems * ncells_elem

        allocate(xi_grid(refine_level), eta_grid(refine_level))
        allocate(N_eval(nbasis), reordered_coords(nbasis, 2))
        allocate(weights_local(nbasis))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, NGRP))
        allocate(Cells(n_sub_elems, 4))

        ! Define interpolation coordinates in reference space [-1, 1]
        do i = 1, refine_level
            xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / real(p, dp)
            eta_grid(i) = xi_grid(i)
        end do

        ! 1. Interpolate geometry and solution to the sub-grid
        gid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*nbasis
            ! Get actual physical coordinates of the high-order nodes
            do i = 1, nbasis
                reordered_coords(i, :) = mesh%nodes(mesh%elems(ee, i), 1:2)
            end do
            weights_local = mesh%weights(mesh%elems(ee, 1:nbasis))
            
            ! Loop through eta (j) then xi (i) to match standard VTK ordering
            do j = 1, refine_level
                do i = 1, refine_level
                    gid = gid + 1
                    call GetArbitraryBasis(FE, xi_grid(i), eta_grid(j), N_eval, weights_local)
                    
                    Xp(gid,1) = dot_product(N_eval, reordered_coords(:, 1))
                    Xp(gid,2) = dot_product(N_eval, reordered_coords(:, 2))
                    Xp(gid,3) = 0.0_dp
                    
                    do g = 1, NGRP
                        Up(gid, g) = dot_product(N_eval, XPETSc(basep+1:basep+nbasis, g))
                    end do
                end do
            end do
        end do

        ! 2. Generate sub-quad connectivity (Linear quads for visualization)
        cid = 0
        do ee = 1, mesh%n_elems
            ! Start of this element's points in the global Xp array
            basep = (ee-1)*npts_elem 
            do j = 1, refine_level - 1
                do i = 1, refine_level - 1
                    ! Node indices in the gid sequence (1-based for now)
                    n00 = basep + (j-1)*refine_level + i
                    n10 = n00 + 1
                    n01 = basep + j*refine_level + i
                    n11 = n01 + 1
                    
                    cid = cid + 1
                    ! Store as counter-clockwise: bottom-left, bottom-right, top-right, top-left
                    Cells(cid, :) = [n00, n10, n11, n01]
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
        
        if (slash_idx > 0) then
            full_path = trim(filename) // "/" // filename(slash_idx+1:) // trim(sem) // trim(adj) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSN%order)) // ".vtk"
        else
            full_path = trim(filename) // "/" // trim(filename) // trim(sem) // trim(adj) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSN%order)) // ".vtk"
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

        write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*5
        do cid = 1, n_sub_elems
            write(unit_v, '(5I10)') 4, Cells(cid,1)-1, Cells(cid,2)-1, &
                                       Cells(cid,3)-1, Cells(cid,4)-1
        end do

        write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
        do i = 1, n_sub_elems
            write(unit_v, '(I2)') 9 ! VTK_QUAD
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
        
        write(unit_v, '(A)') "SCALARS Pins_ID int 1"
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
        deallocate(xi_grid, eta_grid, N_eval, reordered_coords, Xp, Up, Cells)

    end subroutine export_dfem_vtk

end module m_outpVTK_dfem