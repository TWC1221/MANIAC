#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>

module m_outpVTK
  use m_constants
  use m_types
  use m_petsc
  use m_finite_elements

  use petscsys
  use petscvec
  use petscis

  implicit none
  private
  public :: export_vtk_petsc, export_vtk_pcg, derive_case_nametag

  contains

  pure function derive_case_nametag(filename) result(tag)
    character(len=*), intent(in) :: filename
    character(len=:), allocatable :: tag
    integer :: pos

    pos = index(filename, '/', back=.true.)
    if (pos == 0) pos = index(filename, '\', back=.true.)

    tag = filename(pos+1:)
  end function derive_case_nametag

subroutine export_vtk_pcg(filename, FE, mesh, X_PCG, NGRP, refine_level, use_z, printout)
  character(len=*), intent(in) :: filename
  type(t_finite),   intent(in) :: FE
  type(t_mesh),     intent(in) :: mesh
  type(t_vec),      intent(in) :: X_PCG(:)
  integer,          intent(in) :: NGRP, refine_level
  logical,          intent(in) :: use_z
  logical, intent(in), optional :: printout

  integer :: ee, g, i, j, k, a, unit_v, gid, cid, basep, p, nbasis
  integer :: npts_elem, ncells_elem, n_sub_nodes, n_sub_elems
  real(dp), allocatable :: xi_grid(:), eta_grid(:), zeta_grid(:), N_eval(:)
  real(dp), allocatable :: Xp(:,:), Up(:,:)
  real(dp), allocatable :: BasisCache(:,:)
  integer,  allocatable :: Cells(:,:) ! 8 nodes for hex
  real(dp) :: xi, eta, zeta
  integer  :: v(8)
  integer  :: pt

  p         = FE%order
  nbasis    = FE%n_basis
  npts_elem = refine_level**3
  ncells_elem = (refine_level-1)**3
  n_sub_nodes = mesh%n_elems * npts_elem
  n_sub_elems = mesh%n_elems * ncells_elem

  allocate(xi_grid(refine_level), eta_grid(refine_level), zeta_grid(refine_level))
  allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, NGRP), Cells(n_sub_elems, 8))
  allocate(BasisCache(npts_elem, nbasis))

  do i = 1, refine_level
    xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
    eta_grid(i) = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
    zeta_grid(i) = -1.0_dp + 2.0_dp * real(i-1, dp) / real(refine_level-1, dp)
  end do

  ! Precompute Basis functions for the refinement grid
  pt = 0
  do k = 1, refine_level
    zeta = zeta_grid(k)
    do j = 1, refine_level
      eta = eta_grid(j)
      do i = 1, refine_level
        xi = xi_grid(i)
        pt = pt + 1
        call GetArbitraryBasis(FE, xi, eta, zeta, BasisCache(pt, :))
      end do
    end do
  end do

  do ee = 1, mesh%n_elems
    basep = (ee-1)*npts_elem
    do pt = 1, npts_elem
      gid = basep + pt
      Xp(gid,1:3) = 0.0_dp
      Up(gid,:)   = 0.0_dp
      do a = 1, nbasis
        Xp(gid,1) = Xp(gid,1) + BasisCache(pt, a) * mesh%nodes(mesh%elems(ee,a), 1)
        Xp(gid,2) = Xp(gid,2) + BasisCache(pt, a) * mesh%nodes(mesh%elems(ee,a), 2)
        Xp(gid,3) = Xp(gid,3) + BasisCache(pt, a) * mesh%nodes(mesh%elems(ee,a), 3)
        do g = 1, NGRP
           Up(gid,g) = Up(gid,g) + BasisCache(pt, a) * X_PCG(g)%vec(mesh%elems(ee,a))
        end do
      end do
    end do
  end do

  if (use_z) Xp(:,3) = Up(:,1)

  cid = 0
  do ee = 1, mesh%n_elems
      basep = (ee-1)*npts_elem
      do k = 1, refine_level-1
        do j = 1, refine_level-1
          do i = 1, refine_level-1
            v(1) = basep + (k-1)*refine_level*refine_level + (j-1)*refine_level + (i-1)
            v(2) = v(1) + 1
            v(3) = v(1) + refine_level + 1
            v(4) = v(1) + refine_level
            v(5) = v(1) + refine_level*refine_level
            v(6) = v(5) + 1
            v(7) = v(5) + refine_level + 1
            v(8) = v(5) + refine_level
            cid = cid + 1
            Cells(cid,:) = v
          end do
        end do
      end do
  end do

  unit_v = 101
  open(unit=unit_v, file=trim(filename), status='replace', action='write')
  write(unit_v, '(A)') "# vtk DataFile Version 3.0"
  write(unit_v, '(A)') "PCG FEM Solution Output"
  write(unit_v, '(A)') "ASCII"
  write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

  write(unit_v, '(A, I10, A)') "POINTS ", n_sub_nodes, " double"
  do gid = 1, n_sub_nodes
    write(unit_v, '(3F20.12)') Xp(gid,1), Xp(gid,2), Xp(gid,3)
  end do

  write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*9
  do cid = 1, n_sub_elems
    write(unit_v, '(I10,8(1X,I10))') 8, Cells(cid,1), Cells(cid,2), Cells(cid,3), Cells(cid,4), Cells(cid,5), Cells(cid,6), Cells(cid,7), Cells(cid,8)
  end do

  write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
  do cid = 1, n_sub_elems
    write(unit_v, '(I10)') 12 ! VTK_HEXAHEDRON
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

  if (present(printout) .and. (printout .eqv. .true. )) write(*,'(A)') " [  I/O   ] :: Visualization saved to " // trim(filename)

  deallocate(xi_grid, eta_grid, zeta_grid, BasisCache, Xp, Up, Cells)
end subroutine export_vtk_pcg

subroutine export_vtk_petsc(filename, FE, mesh, XPETSc, NGRP, refine_level, use_z, printout)
  use m_petsc 
  implicit none
  character(len=*), intent(in) :: filename
  type(t_finite),   intent(in) :: FE
  type(t_mesh),     intent(in) :: mesh
  Vec,              intent(in) :: XPETSc(:)
  integer,          intent(in) :: NGRP, refine_level
  logical,          intent(in) :: use_z
  logical, intent(in), optional :: printout

  integer :: ee, g, i, j, k, a, unit_v
  integer :: nbasis, npts_elem, ncells_elem
  integer :: n_sub_nodes, n_sub_elems
  integer :: gid, cid, basep, node_idx
  integer :: p, pt

  real(dp), allocatable :: xi_grid(:), eta_grid(:), zeta_grid(:)
  real(dp), allocatable :: BasisCache(:,:)
  real(dp) :: xi, eta, zeta

  real(dp), allocatable :: Xp(:,:), Up(:,:)
  integer,  allocatable :: Cells(:,:)
  PetscScalar, pointer :: p_sol(:)
  PetscErrorCode :: ierr

  integer :: n00, v(8)

  ! --- 1. Setup Dimensions ---
  p         = FE%order
  nbasis    = FE%n_basis
  npts_elem = refine_level**3
  ncells_elem = (refine_level-1)**3
  n_sub_nodes = mesh%n_elems * npts_elem
  n_sub_elems = mesh%n_elems * ncells_elem

  allocate(xi_grid(refine_level), eta_grid(refine_level), zeta_grid(refine_level))
  allocate(BasisCache(npts_elem, nbasis))
  allocate(Xp(n_sub_nodes, 3))
  allocate(Up(n_sub_nodes, NGRP))
  allocate(Cells(n_sub_elems, 8))

  ! --- 2. Create Reference Grid ---
  do i = 1, refine_level
    xi_grid(i)  = -1.0_dp + 2.0_dp * real(i-1, dp) / max(1.0_dp, real(refine_level-1, dp))
    eta_grid(i) = -1.0_dp + 2.0_dp * real(i-1, dp) / max(1.0_dp, real(refine_level-1, dp))
    zeta_grid(i)= -1.0_dp + 2.0_dp * real(i-1, dp) / max(1.0_dp, real(refine_level-1, dp))
  end do

  ! Precompute Basis functions at sub-grid points
  pt = 0
  do k = 1, refine_level
    do j = 1, refine_level
      do i = 1, refine_level
        pt = pt + 1
        call GetArbitraryBasis(FE, xi_grid(i), eta_grid(j), zeta_grid(k), BasisCache(pt, :))
      end do
    end do
  end do

  ! --- 3. Interpolate Geometry ---
  do ee = 1, mesh%n_elems
    basep = (ee-1)*npts_elem
    do pt = 1, npts_elem
      gid = basep + pt
      Xp(gid,1:3) = 0.0_dp
      do a = 1, nbasis
        node_idx = mesh%elems(ee,a)
        Xp(gid,1) = Xp(gid,1) + BasisCache(pt, a) * mesh%nodes(node_idx, 1)
        Xp(gid,2) = Xp(gid,2) + BasisCache(pt, a) * mesh%nodes(node_idx, 2)
        Xp(gid,3) = Xp(gid,3) + BasisCache(pt, a) * mesh%nodes(node_idx, 3)
      end do
    end do
  end do

  ! --- 4. Interpolate Solution ---
  do g = 1, NGRP
    PetscCall(VecGetArrayRead(XPETSc(g), p_sol, ierr))
    do ee = 1, mesh%n_elems
      basep = (ee-1)*npts_elem
      do pt = 1, npts_elem
        gid = basep + pt
        Up(gid,g) = 0.0_dp
        do a = 1, nbasis
          node_idx = mesh%elems(ee,a)
          Up(gid,g) = Up(gid,g) + BasisCache(pt, a) * real(p_sol(node_idx), dp)
        end do
      end do
    end do
    PetscCall(VecRestoreArrayRead(XPETSc(g), p_sol, ierr)) ! Fixed parenthesis
  end do

  ! --- 5. Assemble Sub-Cells Connectivity ---
  cid = 0
  do ee = 1, mesh%n_elems
      basep = (ee-1)*npts_elem
      do k = 1, refine_level-1
        do j = 1, refine_level-1
          do i = 1, refine_level-1
            ! Local indices within the npts_elem block
            ! VTK uses 0-based indexing for node references! 
            ! So we subtract 1 from the absolute node position
            n00 = basep + (k-1)*refine_level*refine_level + (j-1)*refine_level + (i-1)
            
            Cells(cid+1,1) = n00
            Cells(cid+1,2) = n00 + 1
            Cells(cid+1,3) = n00 + refine_level + 1
            Cells(cid+1,4) = n00 + refine_level
            Cells(cid+1,5) = n00 + refine_level*refine_level
            Cells(cid+1,6) = Cells(cid+1,5) + 1
            Cells(cid+1,7) = Cells(cid+1,5) + refine_level + 1
            Cells(cid+1,8) = Cells(cid+1,5) + refine_level
            cid = cid + 1
          end do
        end do
      end do
  end do

  unit_v = 101
  open(unit=unit_v, file=trim(filename), status='replace', action='write')

  write(unit_v, '(A)') "# vtk DataFile Version 3.0"
  write(unit_v, '(A)') "FEM Solution Output"
  write(unit_v, '(A)') "ASCII"
  write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

  write(unit_v, '(A, I10, A)') "POINTS ", n_sub_nodes, " double"
  do gid = 1, n_sub_nodes
    write(unit_v, '(3E16.8)') Xp(gid,1), Xp(gid,2), Xp(gid,3)
  end do

  write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*9
  do cid = 1, n_sub_elems
    write(unit_v, '(9I10)') 8, Cells(cid,1), Cells(cid,2), Cells(cid,3), Cells(cid,4), &
                               Cells(cid,5), Cells(cid,6), Cells(cid,7), Cells(cid,8)
  end do

  write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
  do cid = 1, n_sub_elems
    write(unit_v, '(I2)') 12   ! VTK_HEXAHEDRON
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
      write(unit_v, '(E16.8)') Up(gid,g)
    end do
  end do

  close(unit_v)
  deallocate(xi_grid, eta_grid, zeta_grid, BasisCache, Xp, Up, Cells)

  if (present(printout)) then
    if (printout) write(*,'(A)') " [  I/O   ] :: Visualization saved to " // trim(filename)
  endif

end subroutine export_vtk_petsc

subroutine GetArbitraryBasis(FE, xi, eta, zeta, N)
  type(t_finite), intent(in) :: FE
  real(dp), intent(in) :: xi, eta, zeta
  real(dp), intent(out) :: N(:)
  integer :: i, j, k, nn

  nn = 1
  N = 0.0_dp
  do k = 1, FE%order + 1
    do j = 1, FE%order + 1
      do i = 1, FE%order + 1
        N(FE%p(nn)) = FE_basis_1D(FE, FE%order, i, xi) * FE_basis_1D(FE, FE%order, j, eta) * FE_basis_1D(FE, FE%order, k, zeta)
        nn = nn + 1
      end do
    end do
  end do
end subroutine GetArbitraryBasis

end module m_outpVTK