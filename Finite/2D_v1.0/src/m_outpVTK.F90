#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

module m_outpVTK
  use m_constants
  use m_types
  use m_petsc
  use m_finite_elements

  use petscsys
  use petscvec

  implicit none
  private
  public :: export_vtk, export_vtk_pcg

contains
  
subroutine export_vtk_pcg(filename, FE, mesh, X_PCG, NGRP, refine_level, use_z)
    character(len=*), intent(in) :: filename
    type(t_finite),   intent(in) :: FE
    type(t_mesh),     intent(in) :: mesh
    type(t_vec),      intent(in) :: X_PCG(:) ! Using your custom type
    integer,          intent(in) :: NGRP, refine_level
    logical,          intent(in) :: use_z

    integer :: ee, g, i, j, a, unit_v, gid, cid, basep, p, nbasis
    integer :: npts_elem, ncells_elem, n_sub_nodes, n_sub_elems
    real(dp), allocatable :: xi_grid(:), eta_grid(:), N_eval(:)
    real(dp), allocatable :: Xp(:,:), Up(:,:)
    integer,  allocatable :: Cells(:,:)
    real(dp) :: xi, eta
    integer  :: v(4)
    real(dp) :: x(4), y(4)

    p         = FE%order
    nbasis    = FE%n_basis
    npts_elem = refine_level * refine_level
    ncells_elem = (refine_level-1) * (refine_level-1)
    n_sub_nodes = mesh%n_elems * npts_elem
    n_sub_elems = mesh%n_elems * ncells_elem

    allocate(xi_grid(refine_level), eta_grid(refine_level))
    allocate(N_eval(nbasis), Xp(n_sub_nodes, 3), Up(n_sub_nodes, NGRP), Cells(n_sub_elems, 4))

    do i = 1, refine_level
      xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
      eta_grid(i) = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
    end do

    do ee = 1, mesh%n_elems
      basep = (ee-1)*npts_elem
      do j = 1, refine_level
        eta = eta_grid(j)
        do i = 1, refine_level
          xi = xi_grid(i)
          gid = basep + (j-1)*refine_level + i
          call GetArbitraryBasis(FE, xi, eta, N_eval)

          Xp(gid,1:3) = 0.0_dp
          do a = 1, nbasis
            Xp(gid,1) = Xp(gid,1) + N_eval(a) * mesh%nodes(mesh%elems(ee,a), 1)
            Xp(gid,2) = Xp(gid,2) + N_eval(a) * mesh%nodes(mesh%elems(ee,a), 2)
          end do

          do g = 1, NGRP
            Up(gid,g) = 0.0_dp
            do a = 1, nbasis
              Up(gid,g) = Up(gid,g) + N_eval(a) * X_PCG(g)%vec(mesh%elems(ee,a))
            end do
          end do
        end do
      end do
    end do

    if (use_z) Xp(:,3) = Up(:,1)

    cid = 0
    do ee = 1, mesh%n_elems
        basep = (ee-1)*npts_elem
        do j = 1, refine_level-1
            do i = 1, refine_level-1
                v(1) = basep + (j-1)*refine_level + (i-1)
                v(2) = v(1) + 1
                v(3) = v(1) + refine_level + 1
                v(4) = v(1) + refine_level
                
                x = [Xp(v(1)+1,1), Xp(v(2)+1,1), Xp(v(3)+1,1), Xp(v(4)+1,1)]
                y = [Xp(v(1)+1,2), Xp(v(2)+1,2), Xp(v(3)+1,2), Xp(v(4)+1,2)]
                
                call rotate_start_bottom_left(v, x, y)
                if (signed_area(v, Xp) < 0.0_dp) call enforce_ccw(v)
                
                cid = cid + 1
                Cells(cid,:) = v
            end do
        end do
    end do

    ! 4) Write VTK File
    unit_v = 101
    open(unit=unit_v, file=trim(filename)//".vtk", status='replace', action='write')
    write(unit_v, '(A)') "# vtk DataFile Version 3.0"
    write(unit_v, '(A)') "PCG FEM Solution Output"
    write(unit_v, '(A)') "ASCII"
    write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

    write(unit_v, '(A, I10, A)') "POINTS ", n_sub_nodes, " double"
    do gid = 1, n_sub_nodes
      write(unit_v, '(3F20.12)') Xp(gid,1), Xp(gid,2), Xp(gid,3)
    end do

    write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*5
    do cid = 1, n_sub_elems
      write(unit_v, '(I10,4(1X,I10))') 4, Cells(cid,1), Cells(cid,2), Cells(cid,3), Cells(cid,4)
    end do

    write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
    do cid = 1, n_sub_elems
      write(unit_v, '(I10)') 9
    end do

    write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
    write(unit_v, '(A)') "SCALARS Material_ID int 1"
    write(unit_v, '(A)') "LOOKUP_TABLE default"
    do ee = 1, mesh%n_elems
      do i = 1, ncells_elem
        write(unit_v, '(I10)') mesh%mats(ee)
      end do
    end do

    write(unit_v, '(A, I10)') "POINT_DATA ", n_sub_nodes
    do g = 1, NGRP
      write(unit_v, '(A, I0)') "SCALARS Flux_Group_", g
      write(unit_v, '(A)') "double 1"
      write(unit_v, '(A)') "LOOKUP_TABLE default"
      do gid = 1, n_sub_nodes
        write(unit_v, '(F20.12)') Up(gid,g)
      end do
    end do
    close(unit_v)

    deallocate(xi_grid, eta_grid, N_eval, Xp, Up, Cells)
  end subroutine export_vtk_pcg

  subroutine export_vtk(filename, FE, mesh, XPETSc, NGRP, refine_level, use_z)
    character(len=*), intent(in) :: filename
    type(t_finite),   intent(in) :: FE
    type(t_mesh),     intent(in) :: mesh
    Vec,              intent(in) :: XPETSc(:)
    integer,          intent(in) :: NGRP, refine_level
    logical,          intent(in) :: use_z

    integer :: ee, g, i, j, a, unit_v
    integer :: nbasis, npts_elem, ncells_elem
    integer :: n_sub_nodes, n_sub_elems
    integer :: gid, cid, basep, basec
    integer :: p

    real(dp), allocatable :: xi_grid(:), eta_grid(:)
    real(dp) :: xi, eta

    real(dp), allocatable :: N_eval(:)
    real(dp), allocatable :: Xp(:,:), Up(:,:)  ! points (x,y,z), point data
    integer,  allocatable :: Cells(:,:)        ! 4 point indices per quad
    PetscScalar, pointer :: p_sol(:)
    PetscErrorCode :: ierr

    integer :: n00, n10, n11, n01, v(4), k
    real(dp) :: x(4), y(4)

    p         = FE%order
    nbasis    = FE%n_basis
    npts_elem = refine_level * refine_level
    ncells_elem = (refine_level-1) * (refine_level-1)
    n_sub_nodes = mesh%n_elems * npts_elem
    n_sub_elems = mesh%n_elems * ncells_elem

    allocate(xi_grid(refine_level), eta_grid(refine_level))
    allocate(N_eval(nbasis))
    allocate(Xp(n_sub_nodes, 3))
    allocate(Up(n_sub_nodes, NGRP))
    allocate(Cells(n_sub_elems, 4))

    do i = 1, refine_level
      xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
      eta_grid(i) = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
    end do

    gid = 0
    do ee = 1, mesh%n_elems
      do j = 1, refine_level
        eta = eta_grid(j)
        do i = 1, refine_level
          xi = xi_grid(i)
          gid = gid + 1
          call GetArbitraryBasis(FE, xi, eta, N_eval)

          Xp(gid,1) = 0.0_dp
          Xp(gid,2) = 0.0_dp
          Xp(gid,3) = 0.0_dp
          do a = 1, nbasis
            Xp(gid,1) = Xp(gid,1) + N_eval(a) * mesh%nodes(mesh%elems(ee,a), 1)
            Xp(gid,2) = Xp(gid,2) + N_eval(a) * mesh%nodes(mesh%elems(ee,a), 2)
          end do
        end do
      end do
      do g = 1, NGRP
        PetscCall(VecGetArrayRead(XPETSc(g), p_sol, ierr))
        basep = (ee-1)*npts_elem   ! 0-based base for this element’s point block
        do j = 1, refine_level
          eta = eta_grid(j)
          do i = 1, refine_level
            xi = xi_grid(i)
            call GetArbitraryBasis(FE, xi, eta, N_eval)
            gid = basep + (j-1)*refine_level + i         ! 1-based gid
            Up(gid,g) = 0.0_dp
            do a = 1, nbasis
              Up(gid,g) = Up(gid,g) + N_eval(a) * real(p_sol(mesh%elems(ee,a)), dp)
            end do
          end do
        end do
        PetscCall(VecRestoreArrayRead(XPETSc(g), p_sol, ierr))
      end do
    end do

    ! Optionally lift the surface: z := Up(:,1)
    if (use_z) then
      do gid = 1, n_sub_nodes
        Xp(gid,3) = Up(gid,1)
      end do
    end if

    ! 4) Build cell connectivity with robust rotation to bottom-left & CCW
    cid = 0
    do ee = 1, mesh%n_elems
        basep = (ee-1)*npts_elem
        do j = 1, refine_level-1
            do i = 1, refine_level-1
            ! Zero-based local indices in this element block (VTK uses 0-based)
            n00 = basep + (j-1)*refine_level + (i-1)   ! lower-left in param space
            n10 = n00 + 1                               ! +xi
            n01 = n00 + refine_level                    ! +eta
            n11 = n01 + 1                               ! +xi,+eta

            ! Convert to 1-based for array access into Xp (Xp is 1-based)
            x(1) = Xp(n00+1,1); y(1) = Xp(n00+1,2)
            x(2) = Xp(n10+1,1); y(2) = Xp(n10+1,2)
            x(3) = Xp(n11+1,1); y(3) = Xp(n11+1,2)
            x(4) = Xp(n01+1,1); y(4) = Xp(n01+1,2)

            ! Start with standard param order: (n00, n10, n11, n01)
            v = [n00, n10, n11, n01]

            ! Rotate so that v(1) is physical bottom-left (min y, then min x)
            call rotate_start_bottom_left(v, x, y)

            ! Enforce CCW orientation; if signed area < 0, reverse two vertices
            if (signed_area(v, Xp) < 0.0_dp) then
                call enforce_ccw(v)
                write(*, '(A, I6, A, 2(F8.4, A, F8.4, A))') &
                ">>> Element ", ee, " reoriented to CCW. Nodes: (", &
                x(1), ",", y(1), ") and (", x(4), ",", y(4), ")"
            end if

            cid = cid + 1
            Cells(cid,:) = v
            end do
        end do
    end do

    unit_v = 101
    open(unit=unit_v, file=trim(filename)//".vtk", status='replace', action='write')

    write(unit_v, '(A)') "# vtk DataFile Version 3.0"
    write(unit_v, '(A)') "Interpolated FEM Solution"
    write(unit_v, '(A)') "ASCII"
    write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

    write(unit_v, '(A, I10, A)') "POINTS ", n_sub_nodes, " double"
    do gid = 1, n_sub_nodes
      write(unit_v, '(3F20.12)') Xp(gid,1), Xp(gid,2), Xp(gid,3)
    end do

    write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*(4+1)
    do cid = 1, n_sub_elems
      write(unit_v, '(I10,4(1X,I10))') 4, Cells(cid,1), Cells(cid,2), Cells(cid,3), Cells(cid,4)
    end do

    write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
    do cid = 1, n_sub_elems
      write(unit_v, '(I10)') 9   ! VTK_QUAD
    end do

    ! 5) Cell Data: Material IDs
    write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
    write(unit_v, '(A)') "SCALARS Material_ID int 1"
    write(unit_v, '(A)') "LOOKUP_TABLE default"
    
    do ee = 1, mesh%n_elems
        ! Every sub-cell created from the parent element gets the same material tag
        do i = 1, ncells_elem
            write(unit_v, '(I10)') mesh%mats(ee)
        end do
    end do

    write(unit_v, '(A, I10)') "POINT_DATA ", n_sub_nodes
    do g = 1, NGRP
      write(unit_v, '(A, I0)') "SCALARS Flux_Group_", g
      write(unit_v, '(A)') "double 1"
      write(unit_v, '(A)') "LOOKUP_TABLE default"
      do gid = 1, n_sub_nodes
        write(unit_v, '(F20.12)') Up(gid,g)
      end do
    end do

    close(unit_v)

    ! 6) Cleanup
    deallocate(xi_grid, eta_grid, N_eval, Xp, Up, Cells)

    print *, ">>> Interpolated visualization saved to ", trim(filename)//".vtk"
  end subroutine export_vtk

  ! Compute shape values with your basis and FE%p permutation
  subroutine GetArbitraryBasis(FE, xi, eta, N)
    type(t_finite), intent(in) :: FE
    real(dp), intent(in) :: xi, eta
    real(dp), intent(out) :: N(:)
    integer :: i, j, nn
    nn = 1
    N = 0.0_dp
    do j = 1, FE%order + 1
      do i = 1, FE%order + 1
        N(FE%p(nn)) = FE_basis_1D(FE, FE%order, i, xi) * FE_basis_1D(FE, FE%order, j, eta)
        nn = nn + 1
      end do
    end do
  end subroutine GetArbitraryBasis

  ! Rotate v so that v(1) is the bottom-left (min y, then min x) among the 4 points
  subroutine rotate_start_bottom_left(v, x, y)
    integer, intent(inout) :: v(4)
    real(dp), intent(in)   :: x(4), y(4)
    integer :: idx_bl, k
    idx_bl = 1
    do k = 2, 4
      if (y(k) < y(idx_bl).or. (y(k) == y(idx_bl).and. x(k) < x(idx_bl))) idx_bl = k
    end do
    if (idx_bl == 1) return
    v = [ v(mod(idx_bl-1,4)+1), v(mod(idx_bl  ,4)+1), v(mod(idx_bl+1,4)+1), v(mod(idx_bl+2,4)+1) ]
  end subroutine rotate_start_bottom_left

  ! Signed area (2*area) of quad corner loop in Xp; positive => CCW
  pure real(dp) function signed_area(v, Xp) result(A)
    integer, intent(in) :: v(4)
    real(dp), intent(in) :: Xp(:, :)
    integer :: k, k2, idx1, idx2
    A = 0.0_dp
    do k = 1, 4
      k2 = merge(1, k+1, k==4)
      idx1 = v(k)+1   ! Xp is 1-based; v holds 0-based indices
      idx2 = v(k2)+1
      A = A + ( Xp(idx1,1)*Xp(idx2,2) - Xp(idx2,1)*Xp(idx1,2) )
    end do
    ! A is actually 2*area; sign is what we need
  end function signed_area

  ! Reverse orientation to CCW by swapping two vertices (keep v(1) fixed)
  subroutine enforce_ccw(v)
    integer, intent(inout) :: v(4)
    integer :: tmp
    tmp = v(2); v(2) = v(4); v(4) = tmp
  end subroutine enforce_ccw

end module m_outpVTK