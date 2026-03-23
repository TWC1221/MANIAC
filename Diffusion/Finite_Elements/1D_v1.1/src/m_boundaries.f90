module m_boundaries 
    use m_constants
    implicit none
    contains

    subroutine apply_bcs_csr(val, col_ind, row_ptr, F, n_nodes, bc_left, par_left, bc_right, par_right, geom_m, x1)
        use m_constants

        implicit none
        real(dp), intent(inout) :: val(:), F(:), x1
        integer,  intent(inout) :: col_ind(:)
        integer,  intent(in)    :: row_ptr(:)
        integer,  intent(in)    :: n_nodes, geom_m
        integer,  intent(in)    :: bc_left, bc_right
        real(dp), intent(in)    :: par_left, par_right
        integer :: left_node, right_node

        left_node  = 1
        right_node = n_nodes

        if (bc_left == BC_PERIODIC.or. bc_right == BC_PERIODIC) then
            call apply_periodic_csr(val, col_ind, row_ptr, F, n_nodes)
            return
        end if

        select case (bc_left)
        case (BC_DIRICHLET)
            call apply_dirichlet_csr(val, col_ind, row_ptr, F, left_node, par_left)
        case (BC_NEUMANN)
            call apply_neumann_csr(F, left_node, par_left)
        case (BC_VACUUM)
            call apply_vacuum_csr(val, col_ind, row_ptr, left_node, geom_m, x1)         
        case (BC_ALBEDO)
            call apply_albedo_csr(val, col_ind, row_ptr, F, left_node, par_left)  
        case (BC_NONE)
        end select

        select case (bc_right)
        case (BC_DIRICHLET)
            call apply_dirichlet_csr(val, col_ind, row_ptr, F, right_node, par_right)
        case (BC_NEUMANN)
            call apply_neumann_csr(F, right_node, par_right)
        case (BC_VACUUM)
            call apply_vacuum_csr(val, col_ind, row_ptr, right_node, geom_m, x1)
        case (BC_ALBEDO)
            call apply_albedo_csr(val, col_ind, row_ptr, F, right_node, par_right)
        case (BC_NONE)
        end select
    end subroutine apply_bcs_csr

    subroutine apply_dirichlet_csr(val, col_ind, row_ptr, F, node, value)
        use m_constants
        implicit none

        real(dp), intent(inout) :: val(:)
        integer,  intent(in)    :: col_ind(:), row_ptr(:)
        real(dp), intent(inout) :: F(:)
        integer,  intent(in)    :: node
        real(dp), intent(in)    :: value

        integer :: ii, kk, row_start, row_end

        do ii = 1, size(F)
            do kk = row_ptr(ii), row_ptr(ii+1)-1
                if (col_ind(kk) == node) then
                if (ii /= node) F(ii) = F(ii) - val(kk) * value
                exit
                end if
            end do
        end do

        row_start = row_ptr(node)
        row_end   = row_ptr(node+1) - 1
        do kk = row_start, row_end
            val(kk) = 0.0_dp
        end do

        do ii = 1, size(F)
            if (ii == node) cycle
            do kk = row_ptr(ii), row_ptr(ii+1)-1
                if (col_ind(kk) == node) then
                val(kk) = 0.0_dp
                exit
                end if
            end do
        end do

        do kk = row_ptr(node), row_ptr(node+1)-1
            if (col_ind(kk) == node) then
                val(kk) = 1.0_dp
                exit
            end if
        end do

        F(node) = value

    end subroutine apply_dirichlet_csr

    subroutine apply_neumann_csr(F, node, j_out)
        use m_constants
        implicit none

        real(dp), intent(inout) :: F(:)
        integer,  intent(in)    :: node
        real(dp), intent(in)    :: j_out  ! = -D * dφ/dn
        F(node) = F(node) - j_out
    end subroutine apply_neumann_csr

    subroutine apply_vacuum_csr(val, col_ind, row_ptr, node, geom_m, x1)
        use m_constants
        implicit none

        real(dp), intent(inout) :: val(:), x1
        integer,  intent(in)    :: col_ind(:), row_ptr(:)
        integer,  intent(in)    :: node, geom_m

        integer :: kk

        do kk = row_ptr(node), row_ptr(node+1)-1
            if (col_ind(kk) == node) then
                val(kk) = val(kk) + 0.5_dp * x1**geom_m
                exit
            end if
        end do
    end subroutine apply_vacuum_csr

    subroutine apply_albedo_csr(val, col_ind, row_ptr, F, node, albedo)
        use m_constants
        implicit none
        real(dp), intent(inout) :: val(:), F(:)
        integer,  intent(in)    :: col_ind(:), row_ptr(:), node
        real(dp), intent(in)    :: albedo

        integer :: kk
        real(dp) :: alpha

        alpha = 0.5_dp * (1.0_dp - albedo)

        do kk = row_ptr(node), row_ptr(node+1)-1
            if (col_ind(kk) == node) then
            val(kk) = val(kk) + alpha
            exit
            end if
        end do
    end subroutine apply_albedo_csr

    subroutine apply_periodic_csr(val, col_ind, row_ptr, F, n_nodes)
        use m_constants
        implicit none
        integer,  intent(in)    :: n_nodes
        real(dp), intent(inout) :: val(:), F(:)
        integer,  intent(inout) :: col_ind(:)
        integer,  intent(in)    :: row_ptr(:)

        integer :: i, kk, jj, startk, endk, idx_slave, idx_master, len, master_node, slave_node
        real(dp) :: vN
        logical  :: found_master

        master_node = 1
        slave_node = n_nodes

        do i = 1, n_nodes
            startk = row_ptr(i)
            endk   = row_ptr(i+1) - 1
            if (endk < startk) cycle
            do kk = startk, endk
            if (col_ind(kk) == slave_node) then
                vN = val(kk)
                found_master =.false.
                do jj = startk, endk
                if (col_ind(jj) == master_node) then
                    val(jj) = val(jj) + vN
                    val(kk) = 0.0_dp
                    found_master =.true.
                    exit
                end if
                end do

                if (.not. found_master) then
                col_ind(kk) = master_node
                end if
            end if
            end do
        end do

        startk = row_ptr(slave_node)
        endk   = row_ptr(slave_node+1) - 1

        idx_slave  = -1
        idx_master = -1
        do kk = startk, endk
            if (col_ind(kk) == slave_node)  idx_slave  = kk
            if (col_ind(kk) == master_node) idx_master = kk
        end do

        do kk = startk, endk
            val(kk) = 0.0_dp
        end do

        if (idx_slave < 0) then
            idx_slave = startk
            col_ind(idx_slave) = slave_node
        end if

        if (idx_master < 0) then
            if (len >= 2) then
            if (idx_slave == startk) then
                idx_master = startk + 1
            else
                idx_master = startk
            end if
            else
            idx_master = idx_slave
            end if
            col_ind(idx_master) = master_node
        end if

        val(idx_slave)  =  1.0_dp
        val(idx_master) = -1.0_dp
        F(slave_node)   =  0.0_dp
    end subroutine apply_periodic_csr
end module