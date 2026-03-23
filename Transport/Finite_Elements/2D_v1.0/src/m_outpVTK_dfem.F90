module m_outpVTK_dfem
    use m_constants
    use m_types
    use m_finite_elements
    use m_GMSH, only: derive_case_nametag

    implicit none
    private
    public :: export_dfem_vtk

contains

    subroutine export_dfem_vtk(filename, mesh, FE, XPETSc, NGRP)
        character(len=*), intent(in) :: filename
        type(t_mesh),     intent(in) :: mesh
        type(t_finite),   intent(in) :: FE
        real(dp),         intent(in) :: XPETSc(:,:) ! (n_dof, n_groups)
        integer,          intent(in) :: NGRP

        integer :: ee, g, i, j, unit_v
        integer :: nbasis, npts_elem, ncells_elem
        integer :: n_sub_nodes, n_sub_elems
        integer :: gid, cid, basep, p
        integer :: refine_level

        real(dp), allocatable :: xi_grid(:), eta_grid(:)
        real(dp) :: xi, eta

        real(dp), allocatable :: N_eval(:)
        real(dp), allocatable :: reordered_coords(:,:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)  ! points (x,y,z), point data
        integer,  allocatable :: Cells(:,:)        ! 4 point indices per quad

        integer :: n00, n10, n11, n01, v(4)
        real(dp) :: x(4), y(4), area

        integer :: face, e2, f2, i_idx
        integer, allocatable :: pts1(:), pts2(:)
        real(dp) :: tol, diff, maxdiff_dir, maxdiff_rev
        integer, allocatable :: PointMap(:)
        integer :: n_unique, kk
        integer, allocatable :: counts(:)
        real(dp), allocatable :: Xp_u(:,:), Up_u(:,:)
        real(dp) :: tol2
        real(dp) :: gx, gy, dx, dy
        logical :: found
        integer :: newc1, newc2, newc3, newc4

        p         = FE%order
        refine_level = p + 1
        nbasis    = FE%n_basis
        npts_elem = refine_level * refine_level
        ncells_elem = (refine_level-1) * (refine_level-1)
        n_sub_nodes = mesh%n_elems * npts_elem
        n_sub_elems = mesh%n_elems * ncells_elem

        allocate(xi_grid(refine_level), eta_grid(refine_level))
        allocate(N_eval(nbasis))
        allocate(reordered_coords(nbasis, 2))
        allocate(Xp(n_sub_nodes, 3))
        allocate(Up(n_sub_nodes, NGRP))
        allocate(Cells(n_sub_elems, 4))

        do i = 1, refine_level
          xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
          eta_grid(i) = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
        end do

        gid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*nbasis

            do i = 1, nbasis
                reordered_coords(i, :) = mesh%nodes(mesh%elems(ee, i), :)
            end do
            
            do j = 1, refine_level
                eta = eta_grid(j)
                do i = 1, refine_level
                    xi = xi_grid(i)
                    gid = gid + 1
                    call GetArbitraryBasis(FE, xi, eta, N_eval)

                    Xp(gid,1) = dot_product(N_eval, reordered_coords(:, 1))
                    Xp(gid,2) = dot_product(N_eval, reordered_coords(:, 2))
                    Xp(gid,3) = 0.0_dp

                    do g = 1, NGRP
                        Up(gid, g) = dot_product(N_eval, XPETSc(basep+1:basep+nbasis, g))
                    end do
                end do
            end do
        end do

        cid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*npts_elem
            do j = 1, refine_level-1
                do i = 1, refine_level-1
                n00 = basep + (j-1)*refine_level + (i-1)   
                n10 = n00 + 1                             
                n01 = n00 + refine_level                   
                n11 = n01 + 1                               

                x(1) = Xp(n00+1,1); y(1) = Xp(n00+1,2)
                x(2) = Xp(n10+1,1); y(2) = Xp(n10+1,2)
                x(3) = Xp(n11+1,1); y(3) = Xp(n11+1,2)
                x(4) = Xp(n01+1,1); y(4) = Xp(n01+1,2)

                area = 0.5_dp * ( (x(1)*y(2) - y(1)*x(2)) + &
                                  (x(2)*y(3) - y(2)*x(3)) + &
                                  (x(3)*y(4) - y(3)*x(4)) + &
                                  (x(4)*y(1) - y(4)*x(1)) )

                if (area < 0.0_dp) then
                    v = [n00, n01, n11, n10]
                else
                    v = [n00, n10, n11, n01]
                end if

                cid = cid + 1
                Cells(cid,:) = v
                end do
            end do
        end do

    tol = 1.0e-10_dp
    allocate(pts1(refine_level), pts2(refine_level))

    do ee = 1, mesh%n_elems
            basep = (ee-1)*npts_elem
            do face = 1, mesh%n_faces_per_elem
                e2 = mesh%face_info(1, face, ee)
                if (e2 > 0 .and. ee < e2) then
                    f2 = mesh%face_info(2, face, ee)

                    select case (face)
                    case (1) ! bottom: j = 1, i = 1..refine_level
                        do i_idx = 1, refine_level
                            pts1(i_idx) = basep + (i_idx-1) + 1
                        end do
                    case (3) ! top: j = refine_level
                        do i_idx = 1, refine_level
                            pts1(i_idx) = basep + (refine_level-1)*refine_level + (i_idx-1) + 1
                        end do
                    case (2) ! right: i = refine_level, j = 1..refine_level
                        do i_idx = 1, refine_level
                            pts1(i_idx) = basep + (i_idx-1)*refine_level + (refine_level-1) + 1
                        end do
                    case (4) ! left: i = 1, j = 1..refine_level
                        do i_idx = 1, refine_level
                            pts1(i_idx) = basep + (i_idx-1)*refine_level + 1
                        end do
                    end select

                    basep = (e2-1)*npts_elem
                    select case (f2)
                    case (1)
                        do i_idx = 1, refine_level
                            pts2(i_idx) = basep + (i_idx-1) + 1
                        end do
                    case (3)
                        do i_idx = 1, refine_level
                            pts2(i_idx) = basep + (refine_level-1)*refine_level + (i_idx-1) + 1
                        end do
                    case (2)
                        do i_idx = 1, refine_level
                            pts2(i_idx) = basep + (i_idx-1)*refine_level + (refine_level-1) + 1
                        end do
                    case (4)
                        do i_idx = 1, refine_level
                            pts2(i_idx) = basep + (i_idx-1)*refine_level + 1
                        end do
                    end select

                    maxdiff_dir = 0.0_dp
                    maxdiff_rev = 0.0_dp
                    do i_idx = 1, refine_level
                        diff = sqrt((Xp(pts1(i_idx),1) - Xp(pts2(i_idx),1))**2 + &
                                    (Xp(pts1(i_idx),2) - Xp(pts2(i_idx),2))**2)
                        if (diff > maxdiff_dir) maxdiff_dir = diff

                        diff = sqrt((Xp(pts1(i_idx),1) - Xp(pts2(refine_level-i_idx+1),1))**2 + &
                                    (Xp(pts1(i_idx),2) - Xp(pts2(refine_level-i_idx+1),2))**2)
                        if (diff > maxdiff_rev) maxdiff_rev = diff
                    end do

                    if (min(maxdiff_dir, maxdiff_rev) > tol) then
                        write(*,'(A,I6,A,I6,A,I2,A,F12.6)') 'MISMATCH: elems ', ee, ' vs ', e2, ' face ', face, ' maxdiff: ', min(maxdiff_dir, maxdiff_rev)
                        write(*,'(A)') '  sample coords (ee_face -> neighbor_face):'
                        do i_idx = 1, refine_level
                            write(*,'(A,I3,2F14.8, A, I3,2F14.8)') '   e', ee, Xp(pts1(i_idx),1), Xp(pts1(i_idx),2), '  e', e2, Xp(pts2(i_idx),1), Xp(pts2(i_idx),2)
                        end do
                    end if

                    basep = (ee-1)*npts_elem
                end if
            end do
        end do
        deallocate(pts1, pts2)

        unit_v = 101
        open(unit=unit_v, file=trim(filename), status='replace', action='write')

        write(unit_v, '(A)') "# vtk DataFile Version 3.0"
        write(unit_v, '(A)') "Interpolated FEM Solution"
        write(unit_v, '(A)') "ASCII"
        write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

        allocate(PointMap(n_sub_nodes))
        allocate(Xp_u(n_sub_nodes,3))
        allocate(Up_u(n_sub_nodes,NGRP))
        allocate(counts(n_sub_nodes))
        PointMap = -1
        counts = 0
        n_unique = 0
        tol2 = (1.0e-10_dp)**2

        do gid = 1, n_sub_nodes
            gx = Xp(gid,1); gy = Xp(gid,2)
            found = .false.
            do kk = 1, n_unique
                dx = gx - Xp_u(kk,1)
                dy = gy - Xp_u(kk,2)
                if (dx*dx + dy*dy <= tol2) then
                    PointMap(gid) = kk - 1  ! zero-based for VTK
                    counts(kk) = counts(kk) + 1
                    do g = 1, NGRP
                        Up_u(kk,g) = Up_u(kk,g) + Up(gid,g)
                    end do
                    found = .true.
                    exit
                end if
            end do

            if (.not. found) then
                n_unique = n_unique + 1
                Xp_u(n_unique,1:3) = Xp(gid,1:3)
                do g = 1, NGRP
                    Up_u(n_unique,g) = Up(gid,g)
                end do
                counts(n_unique) = 1
                PointMap(gid) = n_unique - 1
            end if
        end do

        do kk = 1, n_unique
            if (counts(kk) > 1) then
                do g = 1, NGRP
                    Up_u(kk,g) = Up_u(kk,g) / real(counts(kk), dp)
                end do
            end if
        end do

        write(unit_v, '(A, I10, A)') "POINTS ", n_unique, " double"
        do kk = 1, n_unique
          write(unit_v, '(3F20.12)') Xp_u(kk,1), Xp_u(kk,2), Xp_u(kk,3)
        end do

        write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*(4+1)
                do cid = 1, n_sub_elems
                    newc1 = PointMap(Cells(cid,1) + 1)
          newc2 = PointMap(Cells(cid,2) + 1)
          newc3 = PointMap(Cells(cid,3) + 1)
          newc4 = PointMap(Cells(cid,4) + 1)
          write(unit_v, '(I10,4(1X,I10))') 4, newc1, newc2, newc3, newc4
        end do

        write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
        do cid = 1, n_sub_elems
          write(unit_v, '(I10)') 9   ! VTK_QUAD
        end do

        write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
        write(unit_v, '(A)') "SCALARS Material_ID int 1"
        write(unit_v, '(A)') "LOOKUP_TABLE default"
        
        do ee = 1, mesh%n_elems
            do i = 1, ncells_elem
                write(unit_v, '(I10)') mesh%mats(ee)
            end do
        end do

                write(unit_v, '(A, I10)') "POINT_DATA ", n_unique
                do g = 1, NGRP
                    write(unit_v, '(A, I0)') "SCALARS Flux_Group_", g
                    write(unit_v, '(A)') "double 1"
                    write(unit_v, '(A)') "LOOKUP_TABLE default"
                    do kk = 1, n_unique
                        write(unit_v, '(F20.12)') Up_u(kk,g)
                    end do
                end do

                ! Free deduplication buffers
                if (allocated(PointMap)) deallocate(PointMap)
                if (allocated(Xp_u)) deallocate(Xp_u)
                if (allocated(Up_u)) deallocate(Up_u)
                if (allocated(counts)) deallocate(counts)

        close(unit_v)

        deallocate(xi_grid, eta_grid, N_eval, reordered_coords, Xp, Up, Cells)

        print *, ">>> Interpolated visualization saved to ", trim(filename)
    end subroutine export_dfem_vtk

end module m_outpVTK_dfem
