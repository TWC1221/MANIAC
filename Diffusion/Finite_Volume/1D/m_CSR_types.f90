module CSR_types
    implicit none
    private
    public :: CSRMatrix

    type :: CSRMatrix
        real(8), allocatable :: AA(:)
        integer, allocatable :: JA(:)
        integer, allocatable :: IA(:)
        integer :: N
    end type CSRMatrix

end module CSR_types
