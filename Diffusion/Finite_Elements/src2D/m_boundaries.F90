#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

module m_boundaries
    use m_constants
    use m_types
    use petscsys
    use petscvec
    use petscmat
    use petscksp

    implicit none

    integer, parameter :: BC_VACUUM     = 1  
    integer, parameter :: BC_REFLECTIVE = 2  
    integer, parameter :: BC_DIRICHLET  = 3  
    integer, parameter :: BC_ALBEDO     = 4  

    type :: t_bc_config
        integer  :: mat_id 
        integer  :: bc_type  
        real(dp) :: value    
    end type t_bc_config

contains

subroutine InitialiseBoundaries(bc_config, ID1, BC1, VAL1, ID2, BC2, VAL2)
    type(t_bc_config), intent(inout) :: bc_config(2)
    integer, intent(in) :: ID1, ID2, BC1, BC2
    real(dp), intent(in) :: VAL1, VAL2

    bc_config(1)%mat_id  = ID1              
    bc_config(1)%bc_type = BC1   ! Options: BC_ALBEDO, BC_VACUUM, BC_REFLECTIVE, BC_DIRICHLET
    bc_config(1)%value   = VAL1        

    bc_config(2)%mat_id  = ID2              
    bc_config(2)%bc_type = BC2
    bc_config(2)%value   = VAL2
end subroutine

subroutine apply_bcs(mesh, bc_cfg, A_petsc, B_petsc, A_pcg, B_pcg, solver_choice)
    type(t_mesh), intent(in)          :: mesh
    type(t_bc_config), intent(in)     :: bc_cfg
    Mat, optional, intent(inout)      :: A_petsc
    Vec, optional, intent(inout)      :: B_petsc
    type(t_mat), optional, intent(inout) :: A_pcg
    type(t_vec), optional, intent(inout) :: B_pcg
    integer, intent(in)               :: solver_choice

    integer, allocatable :: bdy_nodes(:)
    logical, allocatable :: mask(:)
    integer :: i, n_found
    real(dp) :: effective_alpha

    allocate(mask(mesh%n_nodes)); mask = .false.
    n_found = 0
    
    ! Identify boundary nodes
    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == bc_cfg%mat_id) then
            mask(mesh%edges(i,1)) = .true.
            mask(mesh%edges(i,2)) = .true.
            if (size(mesh%edges, 2) >= 3) then
                if (mesh%edges(i,3) > 0) mask(mesh%edges(i,3)) = .true.
            end if
            n_found = n_found + 1
        end if
    end do

    if (n_found == 0) return

    bdy_nodes = pack([(i, i=1, mesh%n_nodes)], mask)

    select case (bc_cfg%bc_type)
        case (BC_DIRICHLET) 
            if (solver_choice == 2) then 
                call PETSC_Apply_Dirichlet(A_petsc, B_petsc, bdy_nodes, bc_cfg%value)
            else
                call PCG_Apply_Dirichlet(A_pcg, B_pcg, bdy_nodes, bc_cfg%value)
            end if
            
        case (BC_VACUUM, BC_ALBEDO)
            effective_alpha = 0.0_dp
            if (bc_cfg%bc_type == BC_ALBEDO) effective_alpha = bc_cfg%value

            if (solver_choice == 2) then
                call PETSC_apply_Robin(mesh, bc_cfg%mat_id, effective_alpha, A_petsc)
            else
                call PCG_Apply_Robin(mesh, bc_cfg%mat_id, effective_alpha, A_pcg)
            end if

        case (BC_REFLECTIVE)
    end select

    if (allocated(bdy_nodes)) deallocate(bdy_nodes)
end subroutine apply_bcs

! ==============================================================================
! ALBEDO / VACUUM LOGIC
! ==============================================================================

subroutine PETSC_Apply_Robin(mesh, target_id, alpha, A)
    type(t_mesh), intent(in) :: mesh
    integer, intent(in)      :: target_id
    real(dp), intent(in)     :: alpha
    Mat, intent(inout)       :: A
    integer :: i, n1, n2, n3, ierr
    real(dp) :: L, t_corner, t_mid, beta

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
    
    call MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, ierr)

    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == target_id) then
            n1 = mesh%edges(i,1); n2 = mesh%edges(i,2)
            L = sqrt(sum((mesh%nodes(n1,:) - mesh%nodes(n2,:))**2))
            
            if (size(mesh%edges, 2) >= 3 .and. mesh%edges(i,3) > 0) then
                n3 = mesh%edges(i,3)
                t_corner = (beta * L) / 6.0_dp
                t_mid    = (4.0_dp * beta * L) / 6.0_dp
                call MatSetValue(A, n1-1, n1-1, t_corner, ADD_VALUES, ierr)
                call MatSetValue(A, n2-1, n2-1, t_corner, ADD_VALUES, ierr)
                call MatSetValue(A, n3-1, n3-1, t_mid,    ADD_VALUES, ierr)
            else
                t_corner = 0.5_dp * L * beta
                call MatSetValue(A, n1-1, n1-1, t_corner, ADD_VALUES, ierr)
                call MatSetValue(A, n2-1, n2-1, t_corner, ADD_VALUES, ierr)
            end if
        end if
    end do
    
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
end subroutine

subroutine PCG_Apply_Robin(mesh, target_id, alpha, A)
    type(t_mesh), intent(in) :: mesh
    integer, intent(in)      :: target_id
    real(dp), intent(in)     :: alpha
    type(t_mat), intent(inout) :: A
    integer :: i, n1, n2, n3
    real(dp) :: L, t_corner, t_mid, beta

    beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)

    do i = 1, mesh%n_edges
        if (mesh%edge_mats(i) == target_id) then
            n1 = mesh%edges(i,1); n2 = mesh%edges(i,2)
            L = sqrt(sum((mesh%nodes(n1,:) - mesh%nodes(n2,:))**2))
            
            if (size(mesh%edges, 2) >= 3 .and. mesh%edges(i,3) > 0) then
                n3 = mesh%edges(i,3)
                t_corner = (beta * L) / 6.0_dp
                t_mid    = (4.0_dp * beta * L) / 6.0_dp
                call pcg_add_diag(A, n1, t_corner)
                call pcg_add_diag(A, n2, t_corner)
                call pcg_add_diag(A, n3, t_mid)
            else
                t_corner = 0.5_dp * L * beta
                call pcg_add_diag(A, n1, t_corner)
                call pcg_add_diag(A, n2, t_corner)
            end if
        end if
    end do
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

! ==============================================================================
! DIRICHLET LOGIC
! ==============================================================================

subroutine PETSC_Apply_Dirichlet(A, B, nodes, val)
    Mat, intent(inout) :: A
    Vec, intent(inout) :: B
    integer, intent(in) :: nodes(:)
    real(dp), intent(in) :: val
    integer :: i, ierr
    PetscInt, allocatable :: idx(:)
    Vec :: temp_x  

    allocate(idx(size(nodes)))
    idx = int(nodes - 1, kind=kind(idx))
    
    call VecDuplicate(B, temp_x, ierr)
    call VecZeroEntries(temp_x, ierr)
    
    if (abs(val) > 1e-12) then
        do i = 1, size(idx)
            call VecSetValue(temp_x, idx(i), val, INSERT_VALUES, ierr)
        end do
    end if
    call VecAssemblyBegin(temp_x, ierr)
    call VecAssemblyEnd(temp_x, ierr)

    call MatZeroRowsColumns(A, size(idx), idx, 1.0_dp, temp_x, B, ierr)
    
    call VecDestroy(temp_x, ierr)
    deallocate(idx)
end subroutine

subroutine PCG_Apply_Dirichlet(A, B, nodes, val)
    type(t_mat), intent(inout) :: A
    type(t_vec), intent(inout) :: B
    integer, intent(in)        :: nodes(:)
    real(dp), intent(in)       :: val
    integer :: i, j, k, target

    do i = 1, size(nodes)
        target = nodes(i)
        do j = A%row_ptr(target), A%row_ptr(target+1) - 1
            k = A%col(j)
            if (k /= target) then
                B%vec(k) = B%vec(k) - A%val(j) * val
                call zero_column_entry(A, k, target)
            end if
        end do
        B%vec(target) = val
        do j = A%row_ptr(target), A%row_ptr(target+1) - 1
            A%val(j) = merge(1.0_dp, 0.0_dp, A%col(j) == target)
        end do
    end do
end subroutine

subroutine zero_column_entry(A, row, col_to_zero)
    type(t_mat), intent(inout) :: A
    integer, intent(in) :: row, col_to_zero
    integer :: j
    do j = A%row_ptr(row), A%row_ptr(row+1) - 1
        if (A%col(j) == col_to_zero) then
            A%val(j) = 0.0_dp
            return
        end if
    end do
end subroutine

end module m_boundaries