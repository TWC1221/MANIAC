module m_finite_element
  use m_constants
  implicit none
  private
  public :: t_fem_1d
  public :: init_reference_element_1d, generate_mesh_1d
  public :: lagrange_shape_1d, lagrange_shape_derivative_1d
  public :: SigmaA_of_x, SigmaS_of_x, SigmaT_of_x, SigmaR_of_x, Sigmanuf_of_x, chi_of_x, D_of_x, S_of_x, write_phi_continuous

  type :: t_fem_1d
    integer :: order
    real(dp), allocatable :: xi_nodes(:)
  end type t_fem_1d

contains

  !=========================================================
  ! Reference element: equidistant nodes on [-1, 1]
  !=========================================================
  subroutine init_reference_element_1d(fe, order)
    type(t_fem_1d), intent(out) :: fe
    integer,        intent(in)  :: order
    integer :: ii
    real(dp) :: dxi

    fe%order = order
    allocate(fe%xi_nodes(order + 1))

    dxi = 2.0_dp / real(order, dp)
    do ii = 1, order + 1
      fe%xi_nodes(ii) = -1.0_dp + dxi * real(ii - 1, dp)
    end do
  end subroutine init_reference_element_1d

  !=========================================================
  ! Uniform 1D mesh
  !=========================================================
  subroutine generate_mesh_1d(fe, x0, x1, n_elem, nodes, elements)
    type(t_fem_1d),intent(in) :: fe
    real(dp), intent(in)  :: x0, x1
    integer,  intent(in)  :: n_elem
    real(dp), allocatable, intent(out) :: nodes(:)
    integer,  allocatable, intent(out) :: elements(:,:)
    integer :: ii, jj

    allocate(nodes(n_elem*fe%order + 1))
    allocate(elements(n_elem, fe%order + 1))

    do ii = 1, n_elem*fe%order + 1
      nodes(ii) = x0 + (x1 - x0) * real(ii - 1, dp) / real(n_elem*fe%order, dp)
    end do

    do ii = 1, n_elem
      elements(ii, :) = [( (ii-1)*fe%order + jj, jj = 1, fe%order+1 )]
    end do

  end subroutine generate_mesh_1d

  !=========================================================
  ! Lagrange shape function
  !=========================================================
  pure real(dp) function lagrange_shape_1d(fe, ii, xi) result(Ni)
    type(t_fem_1d), intent(in) :: fe
    integer,        intent(in) :: ii
    real(dp),       intent(in) :: xi
    integer :: kk

    Ni = 1.0_dp
    do kk = 1, fe%order + 1
      if (kk /= ii) then
        Ni = Ni * (xi - fe%xi_nodes(kk)) / &
                    (fe%xi_nodes(ii) - fe%xi_nodes(kk))
      end if
    end do
  end function lagrange_shape_1d

  !=========================================================
  ! Derivative of Lagrange shape function
  !=========================================================
  pure real(dp) function lagrange_shape_derivative_1d(fe, ii, xi) result(dNi)
    type(t_fem_1d), intent(in) :: fe
    integer,        intent(in) :: ii
    real(dp),       intent(in) :: xi
    integer :: kk, jj
    real(dp) :: term

    dNi = 0.0_dp
    do kk = 1, fe%order + 1
      if (kk /= ii) then
        term = 1.0_dp / (fe%xi_nodes(ii) - fe%xi_nodes(kk))
        do jj = 1, fe%order + 1
          if (jj /= ii .and. jj /= kk) then
            term = term * (xi - fe%xi_nodes(jj)) / &
                           (fe%xi_nodes(ii) - fe%xi_nodes(jj))
          end if
        end do
        dNi = dNi + term
      end if
    end do
  end function lagrange_shape_derivative_1d

  !=========================================================
  ! CSR Sparsity Pattern
  !=========================================================
  subroutine build_CSR_sparsity_pattern(elements, n_nodes, n_elem, order, nnz, row_ptr, col_ind)
    implicit none

    integer, intent(in)  :: n_nodes, n_elem, order
    integer, intent(in)  :: elements(n_elem, order+1)

    integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)

    integer :: idx, jdx, ee, ii, jj
    integer :: nnz, kk
    integer, allocatable :: marker(:)

    allocate(row_ptr(n_nodes+1), marker(n_nodes))

    kk = 1
    nnz = 0
    marker = 0
    row_ptr(1) = 1

    do idx = 1, n_nodes
      do ee = 1, n_elem
        do ii = 1, order+1
          if (elements(ee,ii) == idx) then
            do jj = 1, order+1
              jdx = elements(ee,jj)
              if (marker(jdx) == 0) then
                marker(jdx) = 1
                nnz = nnz + 1
              end if
            end do
          end if
        end do
      end do
      row_ptr(idx+1) = nnz + 1
      marker = 0
    end do

    allocate(col_ind(nnz))
    marker = 0

    do idx = 1, n_nodes
      do ee = 1, n_elem
        do ii = 1, order+1
          if (elements(ee,ii) == idx) then
            do jj = 1, order+1
              jdx = elements(ee,jj)
              if (marker(jdx) == 0) then
                marker(jdx) = 1
                col_ind(kk) = jdx
                kk = kk + 1
              end if
            end do
          end if
        end do
      end do
      marker = 0
    end do

    deallocate(marker)

  end subroutine build_csr_sparsity_pattern

  !=========================================================
  ! OUTPUT
  !=========================================================
    
  subroutine write_phi_continuous(filename, fe, elements, x_nodes, MAT, nplot)
    use m_constants
    use m_CSR_types
    implicit none

    character(len=*),        intent(in) :: filename
    type(t_fem_1d),          intent(in) :: fe
    integer,                 intent(in) :: elements(:,:)
    real(dp),                intent(in) :: x_nodes(:)
    type(t_CSRMatrix),       intent(in) :: MAT(:)
    integer,                 intent(in) :: nplot

    integer :: ee, ii, kk, unit, nloc, nnode, gg, G
    real(dp) :: xi, dxi, xm, xr, x, phi
    real(dp), allocatable :: ue(:), xe(:)
    logical, allocatable :: is_boundary(:)

    nloc  = fe%order + 1
    nnode = size(x_nodes)
    G     = size(MAT)

    allocate(ue(nloc), xe(nloc))
    allocate(is_boundary(nnode))
    is_boundary =.false.

    ! Mark "boundary" nodes (element endpoints); keeps your original behavior
    do ee = 1, size(elements,1)
      is_boundary(elements(ee,1))   =.true.
      is_boundary(elements(ee,nloc)) =.true.
    end do

    open(newunit=unit, file=filename, status="replace", action="write")

    do gg = 1, G
      ! Safety check
      if (size(MAT(gg)%u) /= nnode) then
        write(*,*) 'write_phi_continuous: size mismatch for MAT(', gg, ')%u.'
        cycle
      end if

      ! Dataset A: boundary nodes for group gg
      write(unit,'(A)') '# group '//trim(adjustl(to_string(gg)))//' dataset: boundary nodes'
      do ii = 1, nnode
        if (is_boundary(ii)) then
            write(unit,'(ES14.6,1X,ES14.6)') x_nodes(ii), MAT(gg)%u(ii)
        end if
      end do

      write(unit,'(A)') ""
      write(unit,'(A)') ""

      ! Dataset B: all nodal points for group gg
      write(unit,'(A)') '# group '//trim(adjustl(to_string(gg)))//' dataset: all nodal points'
      do ii = 1, nnode
        write(unit,'(ES14.6,1X,ES14.6)') x_nodes(ii), MAT(gg)%u(ii)
      end do

      write(unit,'(A)') ""
      write(unit,'(A)') ""

      ! Dataset C: interpolated curve for group gg
      write(unit,'(A)') '# group '//trim(adjustl(to_string(gg)))//' dataset: interpolated curve'
      do ee = 1, size(elements,1)

        do ii = 1, nloc
            xe(ii) = x_nodes(elements(ee,ii))
            ue(ii) = MAT(gg)%u(elements(ee,ii))
        end do

        xm  = 0.5_dp * (xe(1) + xe(nloc))
        xr  = 0.5_dp * (xe(nloc) - xe(1))
        dxi = 2.0_dp / real(nplot-1, dp)

        do kk = 0, nplot-1
          xi  = -1.0_dp + dxi * real(kk,dp)
          x   = xm + xr * xi
          phi = 0.0_dp
          do ii = 1, nloc
              phi = phi + lagrange_shape_1d(fe, ii, xi) * ue(ii)
          end do
          write(unit,'(ES14.6,1X,ES14.6)') x, phi
        end do

      end do

      write(unit,'(A)') ""
      write(unit,'(A)') ""
    end do

    close(unit)
    deallocate(ue, xe, is_boundary)

  end subroutine write_phi_continuous

  pure function to_string(i) result(s)
    integer, intent(in) :: i
    character(len=16)   :: s
    write(s,'(I0)') i
  end function to_string


  !=========================================================
  ! Physics coefficients
  !=========================================================
  pure real(dp) function SigmaT_of_x(x, gg, G) result(St)
    use m_constants  
    implicit none
    real(dp), intent(in) :: x
    integer, intent(in) :: gg, G
    integer :: gg_prime
    St = SigmaA_of_x(x,gg) + sum([(SigmaS_of_x(x,gg,gg_prime),gg_prime=1,G)])
  end function SigmaT_of_x

  pure real(dp) function SigmaR_of_x(x, gg, G) result(Sr)
    use m_constants  
    implicit none
    real(dp), intent(in) :: x
    integer, intent(in) :: gg, G
    Sr = SigmaT_of_x(x,gg, G) - SigmaS_of_x(x,gg,gg)
  end function SigmaR_of_x
  
  pure real(dp) function SigmaA_of_x(x, gg) result(Sa)
    use m_constants  
    implicit none
    real(dp), intent(in) :: x
    real(dp), allocatable :: Sa_vec(:)
    integer, intent(in) :: gg
    Sa_vec = [0.015, 0.04, 0.12]
    Sa = Sa_vec(gg)
  end function SigmaA_of_x

  pure real(dp) function SigmaS_of_x(x, gg, gg_prime) result(Ss)
    use m_constants  
    implicit none
    real(dp), intent(in) :: x
    integer, intent(in) :: gg, gg_prime
    real(dp), allocatable :: Ss_M(:,:)
    Ss_M = transpose(reshape([0.20d0, 0.05d0, 0.10d0, &
                    0.10d0, 0.25d0, 0.07d0, &
                    0.17d0, 0.1d0, 0.30d0], shape=[3,3]))
    Ss = Ss_M(gg,gg_prime) 
  end function SigmaS_of_x

  pure real(dp) function D_of_x(x, gg, G) result(D)
    use m_constants  
    implicit none
    real(dp), intent(in) :: x
    integer, intent(in) :: gg, G
    real(dp), parameter :: TINY = 1.0e-12_dp
    real(dp) :: St

    St = SigmaT_of_x(x, gg, G)
    D  = 1.0_dp / (3.0_dp * max(St, TINY))
  end function D_of_x

  pure real(dp) function S_of_x(x, gg) result(S)
    use m_constants  
    implicit none
    real(dp), intent(in) :: x
    integer, intent(in) :: gg
    S = 1.0_dp
  end function S_of_x

  pure real(dp) function Sigmanuf_of_x(x, gg) result(Sf)
    use m_constants  
    implicit none
    real(dp), intent(in) :: x
    real(dp), allocatable :: Sf_vec(:)
    integer, intent(in) :: gg
    Sf_vec = [0.02, 0.1, 0.35]
    Sf = Sf_vec(gg)
  end function Sigmanuf_of_x

  pure real(dp) function chi_of_x(x, gg) result(c)
    real(dp), intent(in) :: x
    real(dp), allocatable :: chi_vec(:)
    integer,  intent(in) :: gg
    chi_vec = [1.0, 0.0, 0.0]
    c = chi_vec(gg)
  end function chi_of_x

end module m_finite_element
