module m_assembly
  use m_constants
  use m_fem
  use m_materials
  use m_quadrature, only: t_quadrature
  implicit none
  private
  public :: build_group_matrix
  public :: assemble_group_rhs

contains

  subroutine build_group_matrix(mats, elem_mat_map, fem, quad, Nq, dNq, gg, csr_pos, VAL_base, geom)
    type(t_material),  intent(in) :: mats(:)
    integer,           intent(in) :: elem_mat_map(:)
    type(t_fem_1d),    intent(in) :: fem
    type(t_quadrature),intent(in) :: quad
    real(dp),          intent(in) :: Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints)
    integer,           intent(in) :: gg
    integer,           intent(in) :: csr_pos(fem%n_elem, fem%nloc, fem%nloc)
    integer,           intent(in) :: geom
    real(dp),          intent(out):: VAL_base(:)

    integer :: ee, q, ii, jj, pos, mat_id
    real(dp) :: xi, wq, xq, J, geo_w, Dq, Sr
    real(dp), allocatable :: Kel(:,:), xe(:)

    allocate(Kel(fem%nloc,fem%nloc), xe(fem%nloc))
    VAL_base = 0.0_dp

    do ee = 1, fem%n_elem
      mat_id = elem_mat_map(ee)
      Kel = 0.0_dp
      xe(:) = fem%x_nodes(fem%elements(ee,1:fem%nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)
        xq = dot_product(Nq(:,q),  xe)
        J  = dot_product(dNq(:,q), xe)
        geo_w = max(xq, 1.0e-14_dp)**geom
        Dq = mats(mat_id)%D(gg)
        Sr = mats(mat_id)%SigmaR(gg)
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
  end subroutine build_group_matrix

  subroutine assemble_group_rhs(mats, elem_mat_map, fem, quad, Nq, dNq, gg, flux_prev, K_eff_prev, adjoint, F_out, geom)
    type(t_material),  intent(in) :: mats(:)
    integer,           intent(in) :: elem_mat_map(:)
    type(t_fem_1d),    intent(in) :: fem
    type(t_quadrature),intent(in) :: quad
    real(dp),          intent(in) :: Nq(fem%nloc, quad%NoPoints), dNq(fem%nloc, quad%NoPoints)
    integer,           intent(in) :: gg
    real(dp),          intent(in) :: flux_prev(:, :)
    real(dp),          intent(in) :: K_eff_prev
    logical,           intent(in) :: adjoint
    integer,           intent(in) :: geom
    real(dp),          intent(out):: F_out(:)

    integer :: ee, q, ii, jg, mat_id, G
    real(dp) :: xi, wq, xq, J, geo_w, scatter_in_q, fission_prev, fission_coeff
    real(dp), allocatable :: Fel(:), xe(:), ue(:), phi_prev(:)

    G = size(mats(1)%D)
    allocate(Fel(fem%nloc), xe(fem%nloc), ue(fem%nloc), phi_prev(G))
    F_out = 0.0_dp

    do ee = 1, fem%n_elem
      mat_id = elem_mat_map(ee)
      Fel = 0.0_dp
      xe(:) = fem%x_nodes(fem%elements(ee,1:fem%nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)
        xq = dot_product(Nq(:,q),  xe)
        J  = dot_product(dNq(:,q), xe)
        geo_w = max(xq, 1.0e-14_dp)**geom

        do jg = 1, G
          ue = flux_prev(jg, fem%elements(ee,1:fem%nloc))
          phi_prev(jg) = dot_product(ue, Nq(:,q))
        end do

        scatter_in_q = 0.0_dp
        do jg = 1, G
          if (jg /= gg) then
            if (.not. adjoint) then
              scatter_in_q = scatter_in_q + mats(mat_id)%SigmaS(jg, gg) * phi_prev(jg)
            else
              scatter_in_q = scatter_in_q + mats(mat_id)%SigmaS(gg, jg) * phi_prev(jg)
            end if
          end if
        end do

        fission_prev = 0.0_dp
        do jg = 1, G
          if (.not. adjoint) then
            fission_prev = fission_prev + mats(mat_id)%NuSigF(jg) * phi_prev(jg)
          else
            fission_prev = fission_prev + mats(mat_id)%Chi(jg) * phi_prev(jg)
          end if
        end do

        if (.not. adjoint) then
          fission_coeff = mats(mat_id)%Chi(gg)
        else
          fission_coeff = mats(mat_id)%NuSigF(gg)
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
  end subroutine assemble_group_rhs

end module m_assembly