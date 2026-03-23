module m_integrals
  use m_constants
  use m_materials
  use m_fem
  use m_quadrature
  implicit none
  private
  public :: integrate_scalar_over_mesh, l2_error_sq_integrand

contains

  subroutine integrate_scalar_over_mesh(mats, elem_mat_map, fem, quad, Nq, dNq, flux, integrand, result, adjoint, geom)
    type(t_material),  intent(in) :: mats(:)
    type(t_fem_1d),    intent(in) :: fem
    integer,           intent(in) :: elem_mat_map(:)
    type(t_quadrature),intent(in) :: quad
    real(dp),          intent(in) :: Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints)
    real(dp),          intent(in) :: flux(:,:) ! (G, n_nodes)
    interface
        function integrand(mat, xq, phi_vals, adjoint) result(val)
            import t_material, dp
            type(t_material), intent(in) :: mat
            real(dp), intent(in) :: xq
            real(dp), intent(in) :: phi_vals(:)
            logical, intent(in) :: adjoint
            real(dp) :: val
        end function integrand
    end interface
    real(dp),          intent(out):: result
    logical,           intent(in) :: adjoint
    integer,           intent(in) :: geom

    integer :: ee, q, jg, G, mat_id
    real(dp) :: xi, wq, xq, J, geo_w
    real(dp), allocatable :: xe(:), ue(:), phi_vals(:)

    G = size(flux, 1)
    allocate(xe(fem%nloc), ue(fem%nloc), phi_vals(G))
    result = 0.0_dp

    do ee = 1, fem%n_elem
      mat_id = elem_mat_map(ee)
      xe(:) = fem%x_nodes(fem%elements(ee,1:fem%nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)
        xq = dot_product(Nq(:,q),  xe)
        J  = dot_product(dNq(:,q), xe)
        geo_w = max(xq, 1.0e-14_dp)**geom

        do jg = 1, G
          ue = flux(jg, fem%elements(ee,1:fem%nloc))
          phi_vals(jg) = dot_product(ue, Nq(:,q))
        end do
        
        result = result + wq * geo_w * J * integrand(mats(mat_id), xq, phi_vals, adjoint)
      end do
    end do
    deallocate(xe, ue, phi_vals)
  end subroutine integrate_scalar_over_mesh

  function l2_error_sq_integrand(mat, xq, phi_vals, adjoint) result(val)
    type(t_material), intent(in) :: mat
    real(dp), intent(in) :: xq
    real(dp), intent(in) :: phi_vals(:)
    logical, intent(in) :: adjoint
    real(dp) :: val
    real(dp) :: exact

    exact = sin(PI * xq)
    val = (phi_vals(1) - exact)**2
  end function l2_error_sq_integrand

end module m_integrals