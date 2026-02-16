module m_assembly
  use m_constants
  use m_fem
  use m_materials
  use m_quadrature, only: t_quadrature
  implicit none
  private
  public :: build_group_matrix_once
  public :: assemble_group_rhs_only

contains

  ! Build and store the group stiffness matrix VAL_base(gg) once
  subroutine build_group_matrix_once(mats, fem, quad, Nq, dNq, gg, csr_pos, VAL_base)
    type(t_material),  intent(in) :: mats
    type(t_fem_1d),    intent(in) :: fem
    type(t_quadrature),intent(in) :: quad
    real(dp),          intent(in) :: Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints)
    integer,           intent(in) :: gg
    integer,           intent(in) :: csr_pos(fem%n_elem, fem%nloc, fem%nloc)
    real(dp),          intent(out):: VAL_base(:)

    integer :: ee, q, ii, jj, pos
    real(dp) :: xi, wq, xq, J, geo_w, Dq, Sr
    real(dp), allocatable :: Kel(:,:), xe(:)

    allocate(Kel(fem%nloc,fem%nloc), xe(fem%nloc))
    VAL_base = 0.0_dp

    do ee = 1, fem%n_elem
      Kel = 0.0_dp
      xe(:) = fem%x_nodes(fem%elements(ee,1:fem%nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)
        xq = dot_product(Nq(:,q),  xe)
        J  = dot_product(dNq(:,q), xe)
        geo_w = mats%geo_weight(xq)
        Dq = mats%D_of(xq, gg)
        Sr = mats%SigmaR_of(xq, gg)
        do ii = 1, fem%nloc
          do jj = 1, fem%nloc
            Kel(ii,jj) = Kel(ii,jj) + wq * geo_w * ( Dq * (dNq(ii,q)*dNq(jj,q)) / J + Sr * Nq(ii,q) * Nq(jj,q) * J )
          end do
        end do
      end do
      do ii = 1, fem%nloc
        do jj = 1, fem%nloc
          pos = csr_pos(ee,ii,jj)
          VAL_base(pos) = VAL_base(pos) + Kel(ii,jj)
        end do
      end do
    end do

    deallocate(Kel, xe)
  end subroutine build_group_matrix_once

  subroutine assemble_group_rhs_only(mats, fem, quad, Nq, dNq, gg, flux_prev, K_eff_prev, adjoint, F_out)
    type(t_material),  intent(in) :: mats
    type(t_fem_1d),    intent(in) :: fem
    type(t_quadrature),intent(in) :: quad
    real(dp),          intent(in) :: Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints)
    integer,           intent(in) :: gg
    real(dp),          intent(in) :: flux_prev(mats%G, fem%n_nodes)
    real(dp),          intent(in) :: K_eff_prev
    logical,           intent(in) :: adjoint
    real(dp),          intent(out):: F_out(:)

    integer :: ee, q, ii, jg
    real(dp) :: xi, wq, xq, J, geo_w, scatter_in_q, fission_prev, fission_coeff
    real(dp), allocatable :: Fel(:), xe(:), ue(:), phi_prev(:)

    allocate(Fel(fem%nloc), xe(fem%nloc), ue(fem%nloc), phi_prev(mats%G))
    F_out = 0.0_dp

    do ee = 1, fem%n_elem
      Fel = 0.0_dp
      xe(:) = fem%x_nodes(fem%elements(ee,1:fem%nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)
        xq = dot_product(Nq(:,q),  xe)
        J  = dot_product(dNq(:,q), xe)
        geo_w = mats%geo_weight(xq)

        do jg = 1, mats%G
          ue = flux_prev(jg, fem%elements(ee,1:fem%nloc))
          phi_prev(jg) = dot_product(ue, Nq(:,q))
        end do

        scatter_in_q = 0.0_dp
        do jg = 1, mats%G
          if (jg /= gg) then
            if (.not. adjoint) then
              scatter_in_q = scatter_in_q + mats%SigmaS_of(xq, jg, gg) * phi_prev(jg)
            else
              scatter_in_q = scatter_in_q + mats%SigmaS_of(xq, gg, jg) * phi_prev(jg)
            end if
          end if
        end do

        fission_prev = 0.0_dp
        do jg = 1, mats%G
          if (.not. adjoint) then
            fission_prev = fission_prev + mats%SigmanuF_of(xq, jg) * phi_prev(jg)
          else
            fission_prev = fission_prev + mats%chi_of(xq, jg) * phi_prev(jg)
          end if
        end do

        if (.not. adjoint) then
          fission_coeff = mats%chi_of(xq, gg)
        else
          fission_coeff = mats%SigmanuF_of(xq, gg)
        end if

        do ii = 1, fem%nloc
          Fel(ii) = Fel(ii) + wq * geo_w * ( scatter_in_q + fission_coeff * fission_prev / K_eff_prev ) * Nq(ii,q) * J
        end do
      end do

      do ii = 1, fem%nloc
        F_out(fem%elements(ee,ii)) = F_out(fem%elements(ee,ii)) + Fel(ii)
      end do
    end do

    deallocate(Fel, xe, ue, phi_prev)
  end subroutine assemble_group_rhs_only

end module m_assembly