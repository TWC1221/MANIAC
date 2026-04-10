#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

module m_boundaries
    use m_constants
    use m_types
    use m_quadrature
    use m_finite_elements
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
    integer :: i, j, ierr, node_id
    real(dp) :: effective_alpha

    select case (bc_cfg%bc_type)
        case (BC_DIRICHLET) 
            allocate(mask(mesh%n_nodes)); mask = .false.
            do i = 1, mesh%n_edges
                if (mesh%edge_mats(i) == bc_cfg%mat_id) then
                    do j = 1, size(mesh%edges, 2)
                        node_id = mesh%edges(i, j)
                        if (node_id > 0) mask(node_id) = .true.
                    end do
                end if
            end do
            
            if (any(mask)) then
                bdy_nodes = pack([(i, i=1, mesh%n_nodes)], mask)
                if (present(A_petsc)) then 
                    do i = 1, size(bdy_nodes)
                        call MatSetValue(A_petsc, bdy_nodes(i)-1, bdy_nodes(i)-1, PENALTY, ADD_VALUES, ierr)
                    end do
                else
                    do i = 1, size(bdy_nodes)
                        call pcg_add_value(A_pcg, bdy_nodes(i), bdy_nodes(i), PENALTY)
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
    
    integer  :: i, j, k, gp, npe, node_id
    real(dp) :: beta, detJ1D, dx_dxi, dy_dxi
    PetscInt, allocatable :: edge_nodes(:)
    PetscScalar, allocatable :: vals(:)
    PetscErrorCode :: ierr

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
    npe = FE%order + 1 
    
    allocate(edge_nodes(npe), vals(npe*npe))

    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == target_id) then
            do j = 1, npe
                edge_nodes(j) = mesh%edges(i, j) - 1
            end do

            vals = 0.0_dp
            do gp = 1, QuadBound%NoPoints
                dx_dxi = 0.0_dp
                dy_dxi = 0.0_dp
                
                do k = 1, npe
                    node_id = mesh%edges(i, k)
                    if (node_id > 0) then
                        dx_dxi = dx_dxi + mesh%nodes(node_id, 1) * FE%dN_B(gp, k)
                        dy_dxi = dy_dxi + mesh%nodes(node_id, 2) * FE%dN_B(gp, k)
                    end if
                end do
                detJ1D = sqrt(dx_dxi**2 + dy_dxi**2)

                do j = 1, npe
                    do k = 1, npe
                        vals((j-1)*npe + k) = vals((j-1)*npe + k) + &
                            FE%N_B(gp, j) * FE%N_B(gp, k) * &
                            detJ1D * QuadBound%W(gp) * beta
                    end do
                end do
            end do
            
            call MatSetValues(A, npe, edge_nodes, npe, edge_nodes, vals, ADD_VALUES, ierr)
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
    type(t_PCG_CSR),        intent(inout) :: A
    
    integer  :: i, j, k, gp, npe, node_id, node_i, node_j
    real(dp) :: beta, detJ1D, val, dx_dxi, dy_dxi
    real(dp), allocatable :: local_M(:,:)

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
    npe = FE%order + 1 
    
    allocate(local_M(npe, npe))

    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == target_id) then
            local_M = 0.0_dp
            do gp = 1, QuadBound%NoPoints
                dx_dxi = 0.0_dp
                dy_dxi = 0.0_dp
                do k = 1, npe
                    node_id = mesh%edges(i, k)
                    if (node_id > 0) then
                        dx_dxi = dx_dxi + mesh%nodes(node_id, 1) * FE%dN_B(gp, k)
                        dy_dxi = dy_dxi + mesh%nodes(node_id, 2) * FE%dN_B(gp, k)
                    end if
                end do
                detJ1D = sqrt(dx_dxi**2 + dy_dxi**2)
                
                do j = 1, npe
                    do k = 1, npe
                        local_M(j, k) = local_M(j, k) + &
                            FE%N_B(gp, j) * FE%N_B(gp, k) * &
                            detJ1D * QuadBound%W(gp) * beta
                    end do
                end do
            end do
            
            do j = 1, npe
                node_i = mesh%edges(i, j)
                if (node_i <= 0) cycle
                do k = 1, npe
                    node_j = mesh%edges(i, k)
                    if (node_j <= 0) cycle
                    val = local_M(j, k)
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
    integer :: j
    do j = A%row_ptr(row), A%row_ptr(row+1) - 1
        if (A%col(j) == col_idx) then
            A%val(j) = A%val(j) + val
            return
        end if
    end do
    print *, "FATAL ERROR in pcg_add_value: Sparsity pattern mismatch."
    print *, "Attempted to add to non-preallocated entry (", row, ",", col_idx, ")"
    stop
end subroutine

end module m_boundaries