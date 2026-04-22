#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

module m_boundaries
    use m_constants
    use m_types
    use m_quadrature
    use m_basis
    use petscsys
    use petscvec
    use petscmat
    use petscksp

    implicit none

contains

subroutine InitialiseBoundaries(bc_item, ID, BCType, Val)
    type(t_bc_config), intent(out) :: bc_item
    integer,           intent(in)  :: ID, BCType
    real(dp),          intent(in)  :: Val

    bc_item%mat_id  = ID
    bc_item%bc_type = BCType
    bc_item%value   = Val
end subroutine

subroutine apply_bcs(mesh, FE, Quadbound, bc_cfg, A_petsc, A_pcg)
    type(t_mesh), intent(in)          :: mesh
    type(t_finite), intent(in)        :: FE
    type(t_quadrature),intent(in)     :: QuadBound
    type(t_bc_config), intent(in)     :: bc_cfg
    Mat, optional, intent(inout)      :: A_petsc
    type(t_PCG_CSR), optional, intent(inout) :: A_pcg

    integer, allocatable :: bdy_nodes(:)
    logical, allocatable :: mask(:)
    integer :: ii, jj, ierr, node_id
    real(dp) :: effective_alpha

    select case (bc_cfg%bc_type)
        case (BC_DIRICHLET) 
            allocate(mask(mesh%n_nodes)); mask = .false.
            do ii = 1, mesh%n_edges
                if (mesh%edge_mats(ii) == bc_cfg%mat_id) then
                    do jj = 1, size(mesh%edges, 2)
                        node_id = mesh%edges(ii, jj)
                        if (node_id > 0) mask(node_id) = .true.
                    end do
                end if
            end do
            
            if (any(mask)) then
                bdy_nodes = pack([(ii, ii=1, mesh%n_nodes)], mask)
                if (present(A_petsc)) then 
                    do ii = 1, size(bdy_nodes)
                        call MatSetValue(A_petsc, bdy_nodes(ii)-1, bdy_nodes(ii)-1, PENALTY, ADD_VALUES, ierr)
                    end do
                else
                    do ii = 1, size(bdy_nodes)
                        call pcg_add_value(A_pcg, bdy_nodes(ii), bdy_nodes(ii), PENALTY)
                    end do
                end if
                deallocate(bdy_nodes)
            end if
            
        case (BC_VACUUM, BC_ALBEDO)
            effective_alpha = merge(bc_cfg%value, 0.0_dp, bc_cfg%bc_type == BC_ALBEDO)
            if (present(A_petsc)) then
                call PETSC_Apply_Robin(mesh, FE, Quadbound, bc_cfg%mat_id, effective_alpha, A_petsc)
            else
                call PCG_Apply_Robin(mesh, FE, QuadBound, bc_cfg%mat_id, effective_alpha, A_pcg)
            end if
    end select
end subroutine apply_bcs

subroutine PETSC_Apply_Robin(mesh, FE, QuadBound, target_id, alpha, A)    
    type(t_mesh),       intent(in)    :: mesh
    type(t_finite),     intent(in)    :: FE
    type(t_quadrature), intent(in)    :: QuadBound
    integer,            intent(in)    :: target_id
    real(dp),           intent(in)    :: alpha
    Mat,                intent(inout) :: A
    
    integer  :: ii, jj, gp, node_id, ncp, i_span, j_span, row_idx, row_b, col_b
    real(dp) :: beta, detJ2D, u1, u2, v1, v2, xi, eta, dV
    real(dp) :: R(FE%n_basis), dR_dxi(FE%n_basis), dR_deta(FE%n_basis)
    real(dp) :: J(3,2), g_metric(2,2)
    PetscInt, allocatable :: edge_nodes(:)
    PetscScalar, allocatable :: vals(:)
    PetscErrorCode :: ierr

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
    
    allocate(edge_nodes(FE%n_basis), vals(FE%n_basis * FE%n_basis))

    do ii = 1, mesh%n_edges
        if (mesh%edge_mats(ii) == target_id) then
            ncp = mesh%n_cp_edge(ii)
            
            vals(1:ncp*ncp) = 0.0_dp
            do jj = 1, ncp
                edge_nodes(jj) = mesh%edges(ii, jj) - 1
            end do

            do i_span = 1, mesh%n_knots_xi_edge(ii) - 1
                u1 = mesh%edge_knots_xi(ii, i_span); u2 = mesh%edge_knots_xi(ii, i_span+1)
                if (u2 - u1 < 1e-10_dp) cycle
                do j_span = 1, mesh%n_knots_eta_edge(ii) - 1
                    v1 = mesh%edge_knots_eta(ii, j_span); v2 = mesh%edge_knots_eta(ii, j_span+1)
                    if (v2 - v1 < 1e-10_dp) cycle
                    do gp = 1, QuadBound%NoPoints
                        xi  = 0.5_dp * ((u2 - u1) * QuadBound%Xi(gp) + (u2 + u1))
                        eta = 0.5_dp * ((v2 - v1) * QuadBound%Eta(gp) + (v2 + v1))
                        call EvalNURBS2D(FE, ii, mesh, xi, eta, R, dR_dxi, dR_deta)

                        J = 0.0_dp
                        do row_idx = 1, ncp
                            node_id = mesh%edges(ii, row_idx)
                            J(:,1) = J(:,1) + mesh%nodes(node_id, :) * dR_dxi(row_idx)
                            J(:,2) = J(:,2) + mesh%nodes(node_id, :) * dR_deta(row_idx)
                        end do
                        g_metric = matmul(transpose(J), J)
                        detJ2D = sqrt(abs(g_metric(1,1)*g_metric(2,2) - g_metric(1,2)*g_metric(2,1)))
                        dV = detJ2D * QuadBound%W(gp) * beta * 0.25_dp * (u2-u1) * (v2-v1)

                        do row_b = 1, ncp
                            do col_b = 1, ncp
                                vals((row_b-1)*ncp + col_b) = vals((row_b-1)*ncp + col_b) + R(row_b) * R(col_b) * dV
                            end do
                        end do
                    end do
                end do
            end do
            
            call MatSetValues(A, ncp, edge_nodes, ncp, edge_nodes, vals, ADD_VALUES, ierr)
        end if
    end do

    deallocate(edge_nodes, vals)
end subroutine

subroutine PCG_Apply_Robin(mesh, FE, QuadBound, target_id, alpha, A)
    type(t_mesh),       intent(in)    :: mesh
    type(t_finite),     intent(in)    :: FE
    type(t_quadrature), intent(in)    :: QuadBound
    integer,            intent(in)    :: target_id
    real(dp),           intent(in)    :: alpha
    type(t_PCG_CSR),    intent(inout) :: A
    
    integer  :: ii, jj, kk, gp, node_id, node_i, node_j, ncp, i_span, j_span, row_idx
    real(dp) :: beta, detJ2D, val, u1, u2, v1, v2, xi, eta, dV
    real(dp) :: R(FE%n_basis), dR_dxi(FE%n_basis), dR_deta(FE%n_basis)
    real(dp) :: J(3,2), g_metric(2,2)
    real(dp), allocatable :: local_M(:,:)

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
    allocate(local_M(FE%n_basis, FE%n_basis))

    do ii = 1, mesh%n_edges
        if (mesh%edge_mats(ii) == target_id) then
            ncp = mesh%n_cp_edge(ii)
            local_M(1:ncp, 1:ncp) = 0.0_dp

            do i_span = 1, mesh%n_knots_xi_edge(ii) - 1
                u1 = mesh%edge_knots_xi(ii, i_span); u2 = mesh%edge_knots_xi(ii, i_span+1)
                if (u2 - u1 < 1e-10_dp) cycle
                do j_span = 1, mesh%n_knots_eta_edge(ii) - 1
                    v1 = mesh%edge_knots_eta(ii, j_span); v2 = mesh%edge_knots_eta(ii, j_span+1)
                    if (v2 - v1 < 1e-10_dp) cycle
                    do gp = 1, QuadBound%NoPoints
                        xi  = 0.5_dp * ((u2 - u1) * QuadBound%Xi(gp) + (u2 + u1))
                        eta = 0.5_dp * ((v2 - v1) * QuadBound%Eta(gp) + (v2 + v1))
                        call EvalNURBS2D(FE, ii, mesh, xi, eta, R, dR_dxi, dR_deta)

                        J = 0.0_dp
                        do row_idx = 1, ncp
                            node_id = mesh%edges(ii, row_idx)
                            J(:,1) = J(:,1) + mesh%nodes(node_id, :) * dR_dxi(row_idx)
                            J(:,2) = J(:,2) + mesh%nodes(node_id, :) * dR_deta(row_idx)
                        end do
                        g_metric = matmul(transpose(J), J)
                        detJ2D = sqrt(abs(g_metric(1,1)*g_metric(2,2) - g_metric(1,2)*g_metric(2,1)))
                        dV = detJ2D * QuadBound%W(gp) * beta * 0.25_dp * (u2-u1) * (v2-v1)
                        
                        do jj = 1, ncp
                            do kk = 1, ncp
                                local_M(jj, kk) = local_M(jj, kk) + R(jj) * R(kk) * dV
                            end do
                        end do
                    end do
                end do
            end do
            
            do jj = 1, ncp
                node_i = mesh%edges(ii, jj)
                if (node_i <= 0) cycle
                do kk = 1, ncp
                    node_j = mesh%edges(ii, kk)
                    if (node_j <= 0) cycle
                    val = local_M(jj, kk)
                    call pcg_add_value(A, node_i, node_j, val)
                end do
            end do
        end if
    end do

    deallocate(local_M)
end subroutine

subroutine pcg_add_value(A, row, col_idx, val)
    type(t_PCG_CSR), intent(inout) :: A
    integer, intent(in) :: row, col_idx
    real(dp), intent(in) :: val
    integer :: jj
    do jj = A%row_ptr(row), A%row_ptr(row+1) - 1
        if (A%col(jj) == col_idx) then
            A%val(jj) = A%val(jj) + val
            return
        end if
    end do
    print *, "FATAL ERROR in pcg_add_value: Sparsity pattern mismatch."
    print *, "Attempted to add to non-preallocated entry (", row, ",", col_idx, ")"
    stop
end subroutine

end module m_boundaries