module m_output
  use m_constants
  use m_fem,       only: t_fem_1d, lagrange_shape
  use m_CSR_types
  implicit none
  private
  public :: write_phi_continuous

contains

  subroutine write_phi_continuous(filename, fem, MAT, nplot)
    character(len=*), intent(in) :: filename
    type(t_fem_1d),   intent(in) :: fem
    type(t_CSRMatrix),intent(in) :: MAT(:)
    integer,          intent(in) :: nplot

    integer :: ee, ii, kk, unit, nloc, nnode, gg, G
    real(dp) :: xi, dxi, xm, xr, x, phi
    real(dp), allocatable :: ue(:), xe(:)
    logical, allocatable :: is_boundary(:)

    nloc  = fem%nloc
    nnode = fem%n_nodes
    G     = size(MAT)

    allocate(ue(nloc), xe(nloc), is_boundary(nnode))
    is_boundary =.false.
    do ee = 1, fem%n_elem
      is_boundary(fem%elements(ee,1))    =.true.
      is_boundary(fem%elements(ee,nloc)) =.true.
    end do

    open(newunit=unit, file=filename, status="replace", action="write")
    do gg = 1, G
      write(unit,'(A)') '# group '//trim(adjustl(to_string(gg)))//' dataset: boundary nodes'
      do ii = 1, nnode
        if (is_boundary(ii)) write(unit,'(ES14.6,1X,ES14.6)') fem%x_nodes(ii), MAT(gg)%u(ii)
      end do
      write(unit,'(A)') ""; write(unit,'(A)') ""

      write(unit,'(A)') '# group '//trim(adjustl(to_string(gg)))//' dataset: all nodal points'
      do ii = 1, nnode
        write(unit,'(ES14.6,1X,ES14.6)') fem%x_nodes(ii), MAT(gg)%u(ii)
      end do
      write(unit,'(A)') ""; write(unit,'(A)') ""

      write(unit,'(A)') '# group '//trim(adjustl(to_string(gg)))//' dataset: interpolated curve'
      do ee = 1, fem%n_elem
        do ii = 1, nloc
          xe(ii) = fem%x_nodes(fem%elements(ee,ii))
          ue(ii) = MAT(gg)%u(fem%elements(ee,ii))
        end do
        xm  = 0.5_dp * (xe(1) + xe(nloc))
        xr  = 0.5_dp * (xe(nloc) - xe(1))
        dxi = 2.0_dp / real(nplot-1, dp)
        do kk = 0, nplot-1
          xi  = -1.0_dp + dxi * real(kk,dp)
          x   = xm + xr * xi
          phi = dot_product(lagrange_shape(fem, xi), ue)
          write(unit,'(ES14.6,1X,ES14.6)') x, phi
        end do
      end do
      write(unit,'(A)') ""; write(unit,'(A)') ""
    end do
    close(unit)
    deallocate(ue, xe, is_boundary)
  contains
    pure function to_string(i) result(s)
      integer, intent(in) :: i
      character(len=16)   :: s
      write(s,'(I0)') i
    end function to_string
  end subroutine write_phi_continuous

end module m_output