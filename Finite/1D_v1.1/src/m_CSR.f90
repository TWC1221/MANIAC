module m_CSR_types
    use m_constants
    implicit none
    private
    public :: t_CSRMatrix

    type :: t_CSRMatrix
        real(8), allocatable :: val(:), F(:), u(:)
        integer, allocatable :: row_ptr(:)
        integer, allocatable :: col_ind(:)
        integer :: N
    end type t_CSRMatrix
    
    type :: t_quadrature_cache
        integer :: nq, nloc
        real(dp), allocatable :: xi(:), w(:)
        real(dp), allocatable :: Ni(:,:), dNi(:,:)  ! [nloc, nq]
    end type

end module m_CSR_types

module m_csr_utils
  use m_constants
  implicit none
  private
  public :: build_CSR_sparsity_pattern
  public :: build_element_csr_positions

contains

  subroutine build_CSR_sparsity_pattern(elements, n_nodes, n_elem, order, nnz, row_ptr, col_ind)
    integer, intent(in)  :: n_nodes, n_elem, order
    integer, intent(in)  :: elements(n_elem, order+1)
    integer, intent(out) :: nnz
    integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)

    integer :: idx, jdx, ee, ii, jj, kk
    integer, allocatable :: marker(:)

    allocate(row_ptr(n_nodes+1), marker(n_nodes))
    kk = 1; nnz = 0; marker = 0
    row_ptr(1) = 1

    do idx = 1, n_nodes
      do ee = 1, n_elem
        do ii = 1, order+1
          if (elements(ee,ii) == idx) then
            do jj = 1, order+1
              jdx = elements(ee,jj)
              if (marker(jdx) == 0) then
                marker(jdx) = 1
                nnz = nnz + 1
              end if
            end do
          end if
        end do
      end do
      row_ptr(idx+1) = nnz + 1
      marker = 0
    end do

    allocate(col_ind(nnz))
    marker = 0
    kk = 1

    do idx = 1, n_nodes
      do ee = 1, n_elem
        do ii = 1, order+1
          if (elements(ee,ii) == idx) then
            do jj = 1, order+1
              jdx = elements(ee,jj)
              if (marker(jdx) == 0) then
                marker(jdx) = 1
                col_ind(kk) = jdx
                kk = kk + 1
              end if
            end do
          end if
        end do
      end do
      marker = 0
    end do
    deallocate(marker)
  end subroutine build_CSR_sparsity_pattern

  subroutine build_element_csr_positions(n_elem, nloc, elements, row_ptr, col_ind, csr_pos)
    integer, intent(in) :: n_elem, nloc
    integer, intent(in) :: elements(n_elem, nloc)
    integer, intent(in) :: row_ptr(:), col_ind(:)
    integer, intent(out):: csr_pos(n_elem, nloc, nloc)
    integer :: ee, ii, jj, idx, jdx, kk

    do ee = 1, n_elem
      do ii = 1, nloc
        idx = elements(ee,ii)
        do jj = 1, nloc
          jdx = elements(ee,jj)
          csr_pos(ee,ii,jj) = -1
          do kk = row_ptr(idx), row_ptr(idx+1)-1
            if (col_ind(kk) == jdx) then
              csr_pos(ee,ii,jj) = kk
              exit
            end if
          end do
          if (csr_pos(ee,ii,jj) == -1) then
            write(*,*) 'CSR position not found for element ',ee,' (',ii,',',jj,')'
            stop
          end if
        end do
      end do
    end do
  end subroutine build_element_csr_positions

end module m_csr_utils