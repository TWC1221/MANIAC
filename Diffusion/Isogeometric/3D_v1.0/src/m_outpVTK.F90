#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

module m_outpVTK
  use m_constants
  use m_types
  use m_petsc
  use m_basis

  use petscsys
  use petscvec

  implicit none
  private
  public :: derive_case_nametag, export_vtk_petsc, export_vtk_pcg

contains
  
pure function derive_case_nametag(filename) result(tag)
  character(len=*), intent(in) :: filename
  character(len=:), allocatable :: tag
  integer :: pos, dot_pos

  pos = index(filename, '/', back=.true.)
  if (pos == 0) pos = index(filename, '\', back=.true.)

  tag = filename(pos+1:)
  dot_pos = index(tag, '.', back=.true.)
  
  if (dot_pos > 0) then
    tag = tag(1:dot_pos-1) // ".vtk"
  else
    tag = tag // ".vtk"
  end if
end function derive_case_nametag

subroutine export_vtk_pcg(filename, FE, mesh, X_PCG, NGRP, refine_level, use_z)
    character(len=*), intent(in) :: filename
    type(t_finite),   intent(in) :: FE
    type(t_mesh),     intent(in) :: mesh
    type(t_vec),      intent(in) :: X_PCG(:)
    integer,          intent(in) :: NGRP, refine_level
    logical,          intent(in) :: use_z

    integer :: ee, g, ii, jj, kk, a, unit_v, gid, cid, basep, nbasis
    integer :: i_span, j_span, k_span, total_spans
    integer :: n_sub_nodes, n_sub_elems
    real(dp), allocatable :: xi_grid(:), eta_grid(:), zeta_grid(:)
    real(dp) :: R(FE%n_basis), dR_dxi(FE%n_basis), dR_deta(FE%n_basis), dR_dzeta(FE%n_basis)
    real(dp), allocatable :: Xp(:,:), Up(:,:)
    integer,  allocatable :: Cells(:,:)
    real(dp) :: xi, eta, zeta, u1, u2, v1, v2, w1, w2
    integer  :: v(8)

    nbasis    = FE%n_basis

    ! Count total non-zero knot spans across all patches
    total_spans = 0
    do ee = 1, mesh%n_elems
      do ii = 1, size(mesh%knot_vectors_xi, 2) - 1
        if (mesh%knot_vectors_xi(ee, ii+1) - mesh%knot_vectors_xi(ee, ii) > 1e-10_dp) then
          do jj = 1, size(mesh%knot_vectors_eta, 2) - 1
            if (mesh%knot_vectors_eta(ee, jj+1) - mesh%knot_vectors_eta(ee, jj) > 1e-10_dp) then
                do kk = 1, size(mesh%knot_vectors_zeta, 2) - 1
                    if (mesh%knot_vectors_zeta(ee, kk+1) - mesh%knot_vectors_zeta(ee, kk) > 1e-10_dp) then
                        total_spans = total_spans + 1
                    end if
                end do
            end if
          end do
        end if
      end do
    end do

    n_sub_nodes = total_spans * (refine_level ** 3)
    n_sub_elems = total_spans * ((refine_level - 1) ** 3)

    allocate(xi_grid(refine_level), eta_grid(refine_level), zeta_grid(refine_level))
    allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, NGRP), Cells(n_sub_elems, 8))

    do ii = 1, refine_level
      xi_grid(ii)  = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level-1, dp)
      eta_grid(ii) = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level-1, dp)
      zeta_grid(ii) = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level-1, dp)
    end do

    gid = 0
    cid = 0
    do ee = 1, mesh%n_elems
      do i_span = 1, size(mesh%knot_vectors_xi, 2) - 1
        u1 = mesh%knot_vectors_xi(ee, i_span); u2 = mesh%knot_vectors_xi(ee, i_span+1)
        if (u2 - u1 < 1e-10_dp) cycle
        do j_span = 1, size(mesh%knot_vectors_eta, 2) - 1
          v1 = mesh%knot_vectors_eta(ee, j_span); v2 = mesh%knot_vectors_eta(ee, j_span+1)
          if (v2 - v1 < 1e-10_dp) cycle
          do k_span = 1, size(mesh%knot_vectors_zeta, 2) - 1
            w1 = mesh%knot_vectors_zeta(ee, k_span); w2 = mesh%knot_vectors_zeta(ee, k_span+1)
            if (w2 - w1 < 1e-10_dp) cycle

            basep = gid
            do kk = 1, refine_level
              zeta = 0.5_dp * ((w2 - w1) * zeta_grid(kk) + (w2 + w1))
              do jj = 1, refine_level
                eta = 0.5_dp * ((v2 - v1) * eta_grid(jj) + (v2 + v1))
                do ii = 1, refine_level
                  xi = 0.5_dp * ((u2 - u1) * xi_grid(ii) + (u2 + u1))
                  gid = gid + 1
                  call EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R, dR_dxi, dR_deta, dR_dzeta)
                  Xp(gid, 1:3) = 0.0_dp
                  do a = 1, nbasis
                    Xp(gid, 1) = Xp(gid, 1) + R(a) * mesh%nodes(mesh%elems(ee, a), 1)
                    Xp(gid, 2) = Xp(gid, 2) + R(a) * mesh%nodes(mesh%elems(ee, a), 2)
                    Xp(gid, 3) = Xp(gid, 3) + R(a) * mesh%nodes(mesh%elems(ee, a), 3)
                  end do
                  do g = 1, NGRP
                    Up(gid, g) = 0.0_dp
                    do a = 1, nbasis
                      Up(gid, g) = Up(gid, g) + R(a) * X_PCG(g)%vec(mesh%elems(ee, a))
                    end do
                  end do
                end do
              end do
            end do

            do kk = 1, refine_level - 1
              do jj = 1, refine_level - 1
                do ii = 1, refine_level - 1
                  v(1) = basep + (kk-1)*refine_level**2 + (jj-1)*refine_level + (ii-1)
                  v(2) = v(1) + 1
                  v(3) = v(1) + refine_level + 1
                  v(4) = v(1) + refine_level
                  v(5:8) = v(1:4) + refine_level**2
                  cid = cid + 1
                  Cells(cid, :) = v
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    if (use_z) Xp(:,3) = Up(:,1)

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
      write(unit_v, '(I10,8(1X,I10))') 8, (Cells(cid,a), a=1,8)
    end do

    write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
    do cid = 1, n_sub_elems
      write(unit_v, '(I10)') 12
    end do

    write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
    write(unit_v, '(A)') "SCALARS Material_ID int 1"
    write(unit_v, '(A)') "LOOKUP_TABLE default"
    do cid = 1, n_sub_elems
        write(unit_v, '(I10)') 1 ! Simplified matID export for refined grid
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

    print *, ">>> Interpolated visualization saved to ", trim(filename)

    deallocate(xi_grid, eta_grid, zeta_grid, Xp, Up, Cells)
  end subroutine export_vtk_pcg

  subroutine export_vtk_petsc(filename, FE, mesh, XPETSc, NGRP, refine_level, use_z)
    character(len=*), intent(in) :: filename
    type(t_finite),   intent(in) :: FE
    type(t_mesh),     intent(in) :: mesh
    Vec,              intent(in) :: XPETSc(:)
    integer,          intent(in) :: NGRP, refine_level
    logical,          intent(in) :: use_z

    integer :: ee, g, ii, jj, kk, a, unit_v, gid, cid, basep, nbasis
    integer :: i_span, j_span, k_span, total_spans
    integer :: n_sub_nodes, n_sub_elems
    real(dp), allocatable :: xi_grid(:), eta_grid(:), zeta_grid(:)
    real(dp) :: R(FE%n_basis), dR_dxi(FE%n_basis), dR_deta(FE%n_basis), dR_dzeta(FE%n_basis)
    real(dp), allocatable :: Xp(:,:), Up(:,:)
    integer,  allocatable :: Cells(:,:), MatID_p(:)
    PetscScalar, pointer :: p_sol(:)
    PetscErrorCode :: ierr
    real(dp) :: xi, eta, zeta, u1, u2, v1, v2, w1, w2
    integer :: v(8)

    nbasis    = FE%n_basis

    total_spans = 0
    do ee = 1, mesh%n_elems
      do ii = 1, size(mesh%knot_vectors_xi, 2) - 1
        if (mesh%knot_vectors_xi(ee, ii+1) - mesh%knot_vectors_xi(ee, ii) > 1e-10_dp) then
          do jj = 1, size(mesh%knot_vectors_eta, 2) - 1
            if (mesh%knot_vectors_eta(ee, jj+1) - mesh%knot_vectors_eta(ee, jj) > 1e-10_dp) then
              do kk = 1, size(mesh%knot_vectors_zeta, 2) - 1
                if (mesh%knot_vectors_zeta(ee, kk+1) - mesh%knot_vectors_zeta(ee, kk) > 1e-10_dp) then
                  total_spans = total_spans + 1
                end if
              end do
            end if
          end do
        end if
      end do
    end do

    n_sub_nodes = total_spans * (refine_level**3)
    n_sub_elems = total_spans * ((refine_level-1)**3)

    allocate(xi_grid(refine_level), eta_grid(refine_level), zeta_grid(refine_level))
    allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, NGRP), Cells(n_sub_elems, 8), MatID_p(n_sub_elems))

    do ii = 1, refine_level
      xi_grid(ii)  = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level-1, dp)
      eta_grid(ii) = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level-1, dp)
      zeta_grid(ii) = -1.0_dp + 2.0_dp * real(ii-1, dp) / real(refine_level-1, dp)
    end do

    gid = 0
    cid = 0
    do ee = 1, mesh%n_elems
      do i_span = 1, size(mesh%knot_vectors_xi, 2) - 1
        u1 = mesh%knot_vectors_xi(ee, i_span); u2 = mesh%knot_vectors_xi(ee, i_span+1)
        if (u2 - u1 < 1e-10_dp) cycle
        do j_span = 1, size(mesh%knot_vectors_eta, 2) - 1
          v1 = mesh%knot_vectors_eta(ee, j_span); v2 = mesh%knot_vectors_eta(ee, j_span+1)
          if (v2 - v1 < 1e-10_dp) cycle
          do k_span = 1, size(mesh%knot_vectors_zeta, 2) - 1
            w1 = mesh%knot_vectors_zeta(ee, k_span); w2 = mesh%knot_vectors_zeta(ee, k_span+1)
            if (w2 - w1 < 1e-10_dp) cycle

            basep = gid
            do kk = 1, refine_level
              zeta = 0.5_dp * ((w2 - w1) * zeta_grid(kk) + (w2 + w1))
              do jj = 1, refine_level
                eta = 0.5_dp * ((v2 - v1) * eta_grid(jj) + (v2 + v1))
                do ii = 1, refine_level
                  xi = 0.5_dp * ((u2 - u1) * xi_grid(ii) + (u2 + u1))
                  gid = gid + 1
                  call EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R, dR_dxi, dR_deta, dR_dzeta)
                  Xp(gid, 1:3) = 0.0_dp
                  do a = 1, nbasis
                    Xp(gid, 1) = Xp(gid, 1) + R(a) * mesh%nodes(mesh%elems(ee, a), 1)
                    Xp(gid, 2) = Xp(gid, 2) + R(a) * mesh%nodes(mesh%elems(ee, a), 2)
                    Xp(gid, 3) = Xp(gid, 3) + R(a) * mesh%nodes(mesh%elems(ee, a), 3)
                  end do
                  do g = 1, NGRP
                    PetscCall(VecGetArrayRead(XPETSc(g), p_sol, ierr))
                    Up(gid, g) = 0.0_dp
                    do a = 1, nbasis
                      Up(gid, g) = Up(gid, g) + R(a) * real(p_sol(mesh%elems(ee, a)), dp)
                    end do
                    PetscCall(VecRestoreArrayRead(XPETSc(g), p_sol, ierr))
                  end do
                end do
              end do
            end do

            do kk = 1, refine_level - 1
              do jj = 1, refine_level - 1
                do ii = 1, refine_level - 1
                  v(1) = basep + (kk-1)*refine_level**2 + (jj-1)*refine_level + (ii-1)
                  v(2) = v(1) + 1
                  v(3) = v(1) + refine_level + 1
                  v(4) = v(1) + refine_level
                  v(5:8) = v(1:4) + refine_level**2
                  cid = cid + 1
                  Cells(cid, :) = v
                  MatID_p(cid) = mesh%mats(ee)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    if (use_z) Xp(:,3) = Up(:,1)

    unit_v = 101
    open(unit=unit_v, file=trim(filename), status='replace', action='write')
    write(unit_v, '(A)') "# vtk DataFile Version 3.0"
    write(unit_v, '(A)') "PETSc IGA 3D Solution Output"
    write(unit_v, '(A)') "ASCII"
    write(unit_v, '(A)') "DATASET UNSTRUCTURED_GRID"

    write(unit_v, '(A, I10, A)') "POINTS ", n_sub_nodes, " double"
    do gid = 1, n_sub_nodes
      write(unit_v, '(3F20.12)') Xp(gid,1), Xp(gid,2), Xp(gid,3)
    end do

    write(unit_v, '(A, 2I10)') "CELLS ", n_sub_elems, n_sub_elems*9
    do cid = 1, n_sub_elems
      write(unit_v, '(I10,8(1X,I10))') 8, (Cells(cid,a), a=1,8)
    end do

    write(unit_v, '(A, I10)') "CELL_TYPES ", n_sub_elems
    do cid = 1, n_sub_elems
      write(unit_v, '(I10)') 12
    end do

    write(unit_v, '(A, I10)') "CELL_DATA ", n_sub_elems
    write(unit_v, '(A)') "SCALARS Material_ID int 1"
    write(unit_v, '(A)') "LOOKUP_TABLE default"
    do cid = 1, n_sub_elems
        write(unit_v, '(I10)') MatID_p(cid)
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

    deallocate(xi_grid, eta_grid, zeta_grid, Xp, Up, Cells, MatID_p)
    print *, ">>> 3D visualization saved to ", trim(filename)
  end subroutine export_vtk_petsc

end module m_outpVTK
