module m_integrals
  use m_constants
  use m_materials
  use m_fem
  use m_quadrature, only: t_quadrature
  implicit none
  private
  public :: integrate_scalar_over_mesh
  public :: fission_integrand, production_integrand

  abstract interface
    pure function integrand_iface(mats, xq, phi, adjoint) result(val)
      import :: t_material, dp
      type(t_material), intent(in) :: mats
      real(dp),         intent(in) :: xq
      real(dp),         intent(in) :: phi(:)
      logical,          intent(in) :: adjoint
      real(dp) :: val
    end function integrand_iface
  end interface

contains

  subroutine integrate_scalar_over_mesh(mats, fem, quad, Nq, dNq, flux, integrand, total, adjoint)
    type(t_material),  intent(in) :: mats
    type(t_fem_1d),    intent(in) :: fem
    type(t_quadrature),intent(in) :: quad
    real(dp),          intent(in) :: Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints)
    real(dp),          intent(in) :: flux(mats%G, fem%n_nodes)
    procedure(integrand_iface)    :: integrand
    real(dp),          intent(out):: total

    integer :: ee, q, jg
    real(dp) :: xi, wq, xq, J, geo_w, contrib
    logical :: adjoint
    real(dp), allocatable :: xe(:), ue(:), phi(:)

    allocate(xe(fem%nloc), ue(fem%nloc), phi(mats%G))
    total = 0.0_dp

    do ee = 1, fem%n_elem
      xe(:) = fem%x_nodes(fem%elements(ee,1:fem%nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)
        xq = dot_product(Nq(:,q),  xe)
        J  = dot_product(dNq(:,q), xe)
        geo_w = mats%geo_weight(xq)

        do jg = 1, mats%G
          ue = flux(jg, fem%elements(ee,1:fem%nloc))
          phi(jg) = dot_product(ue, Nq(:,q))
        end do

        contrib = integrand(mats, xq, phi, adjoint)
        total   = total + wq * geo_w * contrib * J
      end do
    end do

    deallocate(xe, ue, phi)
  end subroutine integrate_scalar_over_mesh

  pure real(dp) function fission_integrand(mats, xq, phi, adjoint) result(val)
    type(t_material), intent(in) :: mats
    real(dp),         intent(in) :: xq
    real(dp),         intent(in) :: phi(:)
    logical,          intent(in) :: adjoint
    integer :: jg
    val = 0.0_dp
    do jg = 1, size(phi)
      if (.not. adjoint) then
        val = val + mats%SigmanuF_of(xq, jg) * phi(jg)
      else
        val = val + mats%chi_of(xq, jg) * phi(jg)
      end if
    end do
  end function fission_integrand

  pure real(dp) function production_integrand(mats, xq, phi, adjoint) result(val)
    type(t_material), intent(in) :: mats
    real(dp),         intent(in) :: xq
    real(dp),         intent(in) :: phi(:)
    logical,          intent(in) :: adjoint  ! ignored
    integer :: jg
    val = 0.0_dp
    do jg = 1, size(phi)
      val = val + mats%SigmanuF_of(xq, jg) * phi(jg)
    end do
  end function production_integrand

end module m_integrals