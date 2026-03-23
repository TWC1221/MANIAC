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
    fem%order = order
    fem%nloc  = order + 1
    allocate(fem%xi_nodes(fem%nloc))
    if (fem%order > 0) then
      fem%xi_nodes = -1.0_dp + (2.0_dp / real(order, dp)) * real([(ii - 1, ii = 1, fem%nloc)], dp)
    else
      fem%xi_nodes = -1.0_dp
    end if
  end subroutine init_reference_element

  subroutine generate_mesh(fem, x0, x1, n_elem)
    type(t_fem_1d), intent(inout) :: fem
    real(dp),        intent(in)   :: x0, x1
    integer,         intent(in)   :: n_elem
    integer :: ii, jj
    integer :: n_mesh_intervals
    fem%n_elem  = n_elem
    fem%n_nodes = n_elem*fem%order + 1
    allocate(fem%x_nodes(fem%n_nodes))
    allocate(fem%elements(fem%n_elem, fem%nloc))
    n_mesh_intervals = n_elem*fem%order
    if (n_mesh_intervals > 0) then
        fem%x_nodes = x0 + (x1 - x0) * real([(ii - 1, ii = 1, fem%n_nodes)], dp) / real(n_mesh_intervals, dp)
    else
        fem%x_nodes = x0
    end if
    do ii = 1, fem%n_elem
      fem%elements(ii, :) = [( (ii-1)*fem%order + jj, jj = 1, fem%nloc )]
    end do
  end subroutine generate_mesh

  subroutine precompute_shapes(fem, quad, Nq, dNq)
    type(t_fem_1d),    intent(in) :: fem
    type(t_quadrature),intent(in) :: quad
    real(dp), allocatable, intent(out) :: Nq(:,:), dNq(:,:)
    integer :: q
    allocate(Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints))
    do q = 1, quad%NoPoints
      Nq(:,q)  = lagrange_shape(fem, quad%Xi(q))
      dNq(:,q) = lagrange_shape_derivative(fem, quad%Xi(q))
    end do
  end subroutine precompute_shapes

  pure function lagrange_shape(fem, xi) result(Ni)
    type(t_fem_1d), intent(in) :: fem
    real(dp),        intent(in) :: xi
    real(dp) :: Ni(fem%nloc)
    integer :: ii, kk
    do ii = 1, fem%nloc
        Ni(ii) = 1.0_dp
        do kk = 1, fem%nloc
            if (ii == kk) cycle
            Ni(ii) = Ni(ii) * (xi - fem%xi_nodes(kk)) / (fem%xi_nodes(ii) - fem%xi_nodes(kk))
        end do
    end do
  end function lagrange_shape

  pure function lagrange_shape_derivative(fem, xi) result(dNi)
    type(t_fem_1d), intent(in) :: fem
    real(dp),        intent(in) :: xi
    real(dp) :: dNi(fem%nloc)
    integer :: ii, kk, jj
    real(dp) :: term
    do ii = 1, fem%nloc
        dNi(ii) = 0.0_dp
        do kk = 1, fem%nloc
            if (kk == ii) cycle
            term = 1.0_dp / (fem%xi_nodes(ii) - fem%xi_nodes(kk))
            do jj = 1, fem%nloc
                if (jj == ii .or. jj == kk) cycle
                term = term * (xi - fem%xi_nodes(jj)) / (fem%xi_nodes(ii) - fem%xi_nodes(jj))
            end do
            dNi(ii) = dNi(ii) + term
        end do
    end do
  end function lagrange_shape_derivative

end module m_fem