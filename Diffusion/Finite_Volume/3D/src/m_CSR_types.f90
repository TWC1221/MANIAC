module CSR_types
    implicit none
    private
    public :: CSRMatrix, MatMatrix, read_material_file

    type :: CSRMatrix
        real(8), allocatable :: AA(:)
        integer, allocatable :: JA(:)
        integer, allocatable :: IA(:)
        integer :: N
    end type CSRMatrix

    type :: MatMatrix
        integer :: G
        real(8), allocatable :: chi(:)
        real(8), allocatable :: Sigma_s(:,:)
        real(8), allocatable :: nu_Sigma_f(:)
        real(8), allocatable :: Sigma_a(:)
        real(8), allocatable :: EVAL(:)
        real(8), allocatable :: Beta(:)
        real(8), allocatable :: Delayed_Chi(:)
    end type MatMatrix

    contains
    subroutine read_material_file(filename, MATs)
        implicit none
        character(len=*), intent(in) :: filename
        type(MatMatrix), intent(out) :: MATs

        integer :: i, ios, unit, pos
        character(len=256) :: line
        real(8), allocatable :: temp_row(:)
        logical :: eval

        eval = .false.

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print*, 'Error opening material file!'
            stop
        end if

        ! Read G
        do
            read(unit,'(A)', iostat=ios) line
            if (ios /= 0) then
                print*, 'Unexpected end of file while reading G!'
                stop
            end if
            line = adjustl(line)
            if (line(1:5) == '#EVAL') eval = .true.
            if (len_trim(line) == 0) cycle         ! skip empty lines
            if (line(1:1) == '#') cycle           ! skip comments
            ! remove anything before '=' if exists
            pos = index(line,'=')
            if (pos > 0) line = adjustl(line(pos+1:))
            read(line,*) MATs%G
            exit
        end do

        ! Allocate arrays
        allocate(MATs%chi(MATs%G))
        allocate(MATs%Sigma_a(MATs%G))
        allocate(MATs%nu_Sigma_f(MATs%G))
        allocate(MATs%Sigma_s(MATs%G, MATs%G))
        allocate(MATs%EVAL(5))
        allocate(MATs%Beta(6))
        allocate(MATs%Delayed_Chi(6))

        ! Helper to read a vector

        ! Read chi, Sigma_a, nu_Sigma_f
        call read_vector(MATs%chi)
        call read_vector(MATs%Sigma_a)
        call read_vector(MATs%nu_Sigma_f)

        ! Read Sigma_s
        do i = 1, MATs%G
            do
                read(unit,'(A)', iostat=ios) line
                if (ios /= 0) then
                    print*, 'Unexpected EOF reading Sigma_s!'
                    stop
                end if
                line = adjustl(line)
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#') cycle
                pos = index(line,'=')
                if (pos > 0) line = adjustl(line(pos+1:))
                allocate(temp_row(MATs%G))
                read(line,*) temp_row
                MATs%Sigma_s(i,:) = temp_row
                deallocate(temp_row)
                exit
            end do
        end do

        ! call read_vector(MATs%Beta)
        ! call read_vector(MATs%Delayed_Chi)
        if (eval) call read_vector(MATs%EVAL)

        close(unit)
        contains
            subroutine read_vector(arr)
                real(8), intent(out) :: arr(:)
                do
                    read(unit,'(A)', iostat=ios) line
                    if (line(1:5) == '#EVAL') eval = .true.
                    if (ios /= 0) then
                        print*, 'Unexpected EOF reading vector!'
                        stop
                    end if
                    line = adjustl(line)
                    if (len_trim(line) == 0) cycle
                    if (line(1:1) == '#') cycle
                    pos = index(line,'=')
                    if (pos > 0) line = adjustl(line(pos+1:))
                    read(line,*) arr
                    exit
                end do
            end subroutine read_vector
    end subroutine read_material_file

end module CSR_types
