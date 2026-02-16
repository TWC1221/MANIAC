module m_fem
  use m_constants
  use m_quadrature, only: t_quadrature
  implicit none
  private
  public :: t_fem_1d
  public :: init_reference_element, generate_mesh
  public :: precompute_shapes
  public :: lagrange_shape, lagrange_shape_derivative

  type :: t_fem_1d
    integer :: order    = 0
    integer :: n_elem   = 0
    integer :: n_nodes  = 0
    integer :: nloc     = 0
    real(dp), allocatable :: xi_nodes(:)      ! reference nodes [-1,1]
    real(dp), allocatable :: x_nodes(:)       ! mesh coordinates
    integer,  allocatable :: elements(:,:)    ! connectivity [n_elem, nloc]
  end type t_fem_1d

contains

  subroutine init_reference_element(fem, order)
    type(t_fem_1d), intent(inout) :: fem
    integer,        intent(in)    :: order
    integer :: ii
    real(dp) :: dxi
    fem%order = order
    fem%nloc  = order + 1
    allocate(fem%xi_nodes(fem%nloc))
    dxi = 2.0_dp / real(order, dp)
    do ii = 1, fem%nloc
      fem%xi_nodes(ii) = -1.0_dp + dxi * real(ii - 1, dp)
    end do
  end subroutine init_reference_element

  subroutine generate_mesh(fem, x0, x1, n_elem)
    type(t_fem_1d), intent(inout) :: fem
    real(dp),        intent(in)   :: x0, x1
    integer,         intent(in)   :: n_elem
    integer :: ii, jj
    fem%n_elem  = n_elem
    fem%n_nodes = n_elem*fem%order + 1
    allocate(fem%x_nodes(fem%n_nodes))
    allocate(fem%elements(fem%n_elem, fem%nloc))
    do ii = 1, fem%n_nodes
      fem%x_nodes(ii) = x0 + (x1 - x0) * real(ii - 1, dp) / real(n_elem*fem%order, dp)
    end do
    do ii = 1, fem%n_elem
      fem%elements(ii, :) = [( (ii-1)*fem%order + jj, jj = 1, fem%nloc )]
    end do
  end subroutine generate_mesh

  subroutine precompute_shapes(fem, quad, Nq, dNq)
    type(t_fem_1d),    intent(in) :: fem
    type(t_quadrature),intent(in) :: quad
    real(dp), allocatable, intent(out) :: Nq(:,:), dNq(:,:)
    integer :: q, ii
    allocate(Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints))
    do q = 1, quad%NoPoints
      do ii = 1, fem%nloc
        Nq(ii,q)  = lagrange_shape(fem, ii, quad%Xi(q))
        dNq(ii,q) = lagrange_shape_derivative(fem, ii, quad%Xi(q))
      end do
    end do
  end subroutine precompute_shapes

  pure real(dp) function lagrange_shape(fem, ii, xi) result(Ni)
    type(t_fem_1d), intent(in) :: fem
    integer,         intent(in) :: ii
    real(dp),        intent(in) :: xi
    integer :: kk
    Ni = 1.0_dp
    do kk = 1, fem%nloc
      if (kk /= ii) then
        Ni = Ni * (xi - fem%xi_nodes(kk)) / (fem%xi_nodes(ii) - fem%xi_nodes(kk))
      end if
    end do
  end function lagrange_shape

  pure real(dp) function lagrange_shape_derivative(fem, ii, xi) result(dNi)
    type(t_fem_1d), intent(in) :: fem
    integer,         intent(in) :: ii
    real(dp),        intent(in) :: xi
    integer :: kk, jj
    real(dp) :: term
    dNi = 0.0_dp
    do kk = 1, fem%nloc
      if (kk /= ii) then
        term = 1.0_dp / (fem%xi_nodes(ii) - fem%xi_nodes(kk))
        do jj = 1, fem%nloc
          if (jj /= ii.and. jj /= kk) then
            term = term * (xi - fem%xi_nodes(jj)) / (fem%xi_nodes(ii) - fem%xi_nodes(jj))
          end if
        end do
        dNi = dNi + term
      end if
    end do
  end function lagrange_shape_derivative

end module m_fem