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

end module m_CSR_types
