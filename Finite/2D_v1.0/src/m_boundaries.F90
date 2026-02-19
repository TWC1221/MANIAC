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

    integer, parameter :: BC_VACUUM     = 1  
    integer, parameter :: BC_REFLECTIVE = 2  
    integer, parameter :: BC_DIRICHLET  = 3  
    integer, parameter :: BC_ALBEDO     = 4  

    real(dp), parameter :: PENALTY = 1.0e10_dp

    type :: t_bc_config
        integer  :: mat_id 
        integer  :: bc_type  
        real(dp) :: value    
    end type t_bc_config

contains

subroutine InitialiseBoundaries(bc_item, ID, BCType, Val)
    type(t_bc_config), intent(out) :: bc_item
    integer,           intent(in)  :: ID, BCType
    real(dp),          intent(in)  :: Val

    bc_item%mat_id  = ID
    bc_item%bc_type = BCType
    bc_item%value   = Val
end subroutine

subroutine apply_bcs(mesh, FE, Quadbound, bc_cfg, A_petsc, A_pcg, solver_choice)
    type(t_mesh), intent(in)          :: mesh
    type(t_finite), intent(in)        :: FE
    type(t_quadrature),intent(in)     :: QuadBound
    type(t_bc_config), intent(in)     :: bc_cfg
    Mat, optional, intent(inout)      :: A_petsc
    type(t_mat), optional, intent(inout) :: A_pcg
    integer, intent(in)               :: solver_choice

    integer, allocatable :: bdy_nodes(:)
    logical, allocatable :: mask(:)
    integer :: i, j, ierr, node_id
    real(dp) :: effective_alpha

    allocate(mask(mesh%n_nodes)); mask = .false.
    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == bc_cfg%mat_id) then
            do j = 1, size(mesh%edges, 2)
                node_id = mesh%edges(i, j)
                if (node_id > 0) mask(node_id) = .true.
            end do
        end if
    end do
    if (.not. any(mask)) return
    bdy_nodes = pack([(i, i=1, mesh%n_nodes)], mask)

    select case (bc_cfg%bc_type)
        case (BC_DIRICHLET) 
            if (solver_choice == 2) then 
                do i = 1, size(bdy_nodes)
                    call MatSetValue(A_petsc, bdy_nodes(i)-1, bdy_nodes(i)-1, PENALTY, ADD_VALUES, ierr)
                end do
            else
                do i = 1, size(bdy_nodes)
                    call pcg_add_diag(A_pcg, bdy_nodes(i), PENALTY)
                end do
            end if
            
        case (BC_VACUUM, BC_ALBEDO)
            effective_alpha = merge(bc_cfg%value, 0.0_dp, bc_cfg%bc_type == BC_ALBEDO)
            if (solver_choice == 2) then
                call PETSC_Apply_Robin(mesh, FE, Quadbound, bc_cfg%mat_id, effective_alpha, A_petsc)
            else
                call PCG_Apply_Robin(mesh, FE, QuadBound, bc_cfg%mat_id, effective_alpha, A_pcg)
            end if
    end select

    if (allocated(bdy_nodes)) deallocate(bdy_nodes)
end subroutine apply_bcs

subroutine PETSC_Apply_Robin(mesh, FE, QuadBound, target_id, alpha, A)    
    type(t_mesh),       intent(in)    :: mesh
    type(t_finite),     intent(in)    :: FE
    type(t_quadrature), intent(in)    :: QuadBound
    integer,            intent(in)    :: target_id
    real(dp),           intent(in)    :: alpha
    Mat,                intent(inout) :: A
    
    integer  :: i, j, gp, n1, n2, node_id, ierr, npe
    real(dp) :: L, beta, detJ1D
    real(dp), allocatable :: weights(:)
    integer,  allocatable :: edge_map(:)

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
    npe = FE%order + 1 
    
    allocate(weights(npe), edge_map(npe))
    weights = 0.0_dp

    do j = 1, npe
        do gp = 1, QuadBound%NoPoints
            weights(j) = weights(j) + FE%N_B(gp, j) * QuadBound%W(gp)
        end do
    end do

    edge_map(1) = 1
    if (npe > 1) edge_map(2) = npe
    if (npe > 2) then
        do j = 3, npe
            edge_map(j) = j - 1
        end do
    end if

    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == target_id) then
            n1 = mesh%edges(i,1); n2 = mesh%edges(i,2)
            L = sqrt(sum((mesh%nodes(n1,:) - mesh%nodes(n2,:))**2))
            detJ1D = L * 0.5_dp

            do j = 1, npe
                node_id = mesh%edges(i, j)
                if (node_id > 0) then
                    call MatSetValue(A, node_id-1, node_id-1, &
                                     weights(edge_map(j)) * detJ1D * beta, ADD_VALUES, ierr)
                end if
            end do
        end if
    end do

    deallocate(weights, edge_map)
end subroutine

subroutine PCG_Apply_Robin(mesh, FE, QuadBound, target_id, alpha, A)
    type(t_mesh),       intent(in)    :: mesh
    type(t_finite),     intent(in)    :: FE
    type(t_quadrature), intent(in)    :: QuadBound
    integer,            intent(in)    :: target_id
    real(dp),           intent(in)    :: alpha
    type(t_mat),        intent(inout) :: A
    
    integer  :: i, j, gp, n1, n2, node_id, npe
    real(dp) :: L, beta, detJ1D
    real(dp), allocatable :: weights(:)
    integer,  allocatable :: edge_map(:)

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
    npe = FE%order + 1 
    
    allocate(weights(npe), edge_map(npe))
    weights = 0.0_dp

    do j = 1, npe
        do gp = 1, QuadBound%NoPoints
            weights(j) = weights(j) + FE%N_B(gp, j) * QuadBound%W(gp)
        end do
    end do

    edge_map(1) = 1
    if (npe > 1) edge_map(2) = npe
    if (npe > 2) then
        do j = 3, npe
            edge_map(j) = j - 1
        end do
    end if

    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == target_id) then
            n1 = mesh%edges(i,1)
            n2 = mesh%edges(i,2)
            L = sqrt(sum((mesh%nodes(n1,:) - mesh%nodes(n2,:))**2))
            
            detJ1D = L * 0.5_dp

            do j = 1, npe
                node_id = mesh%edges(i, j)
                if (node_id > 0) then
                    call pcg_add_diag(A, node_id, weights(edge_map(j)) * detJ1D * beta)
                end if
            end do
        end if
    end do

    deallocate(weights, edge_map)
end subroutine

subroutine pcg_add_diag(A, node, val)
    type(t_mat), intent(inout) :: A
    integer, intent(in) :: node
    real(dp), intent(in) :: val
    integer :: j
    do j = A%row_ptr(node), A%row_ptr(node+1) - 1
        if (A%col(j) == node) then
            A%val(j) = A%val(j) + val
            return
        end if
    end do
end subroutine

end module m_boundaries