program VTK_NthOrder_Boundary
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    
    type t_edge
        integer :: n(2)
        integer :: id
    end type t_edge

    real(dp), allocatable :: pts(:,:)
    integer, allocatable  :: corner_nodes(:,:)
    integer, allocatable  :: volume_ids(:)
    type(t_edge), allocatable :: all_e(:), b_e(:)
    real(dp) :: xmin, xmax, ymin, ymax, mx, my
    integer  :: npts, ncells, n_b, i, k, csize, n1, n2
    character(len=100) :: infile, outfile
    logical :: exists

    infile = "../input/C5G7vol_3_2_f.vtk"
    outfile = infile

    ! 1. PARSE VTK
    inquire(file=trim(infile), exist=exists)
    if (.not. exists) stop "ERROR: File not found"

    open(10, file=trim(infile), status='old', action='read')
    call parse_points(10, pts, npts)
    call parse_cells(10, corner_nodes, ncells, csize)
    call parse_cell_entity_ids(10, volume_ids, ncells) 
    close(10)

    xmin = minval(pts(:,1)); xmax = maxval(pts(:,1))
    ymin = minval(pts(:,2)); ymax = maxval(pts(:,2))
    
    ! 2. GENERATE EDGES (With safety buffer)
    allocate(all_e((ncells * 4) + 10))
    k = 0
    do i = 1, ncells
        k = k + 1; all_e(k)%n = sort_pair(corner_nodes(i,1), corner_nodes(i,2))
        k = k + 1; all_e(k)%n = sort_pair(corner_nodes(i,2), corner_nodes(i,3))
        k = k + 1; all_e(k)%n = sort_pair(corner_nodes(i,3), corner_nodes(i,4))
        k = k + 1; all_e(k)%n = sort_pair(corner_nodes(i,4), corner_nodes(i,1))
    end do

    call quicksort_edges(all_e, 1, k)

    allocate(b_e(k))
    n_b = 0; i = 1
    do while (i <= k)
        if (i < k .and. all(all_e(i)%n == all_e(i+1)%n)) then
            i = i + 2 
        else
            n_b = n_b + 1
            b_e(n_b) = all_e(i)
            n1 = b_e(n_b)%n(1); n2 = b_e(n_b)%n(2)
            mx = (pts(n1,1) + pts(n2,1)) / 2.0_dp
            my = (pts(n1,2) + pts(n2,2)) / 2.0_dp
            b_e(n_b)%id = label_bc(mx, my, xmin, xmax, ymin, ymax)
            i = i + 1
        end if
    end do

    call write_vtk(outfile, pts, corner_nodes, volume_ids, b_e(1:n_b))
    print *, "DONE: Boundary IDs (101-104) merged with CellEntityIDS."
    
contains

    subroutine parse_cell_entity_ids(u, ids, n)
        integer, intent(in) :: u, n
        integer, allocatable, intent(out) :: ids(:)
        character(len=256) :: line
        integer :: i
        rewind(u)
        allocate(ids(n)); ids = 1
        do while (.true.)
            read(u, '(A)', end=400) line
            if (index(line, "SCALARS CellEntityIds int 1") > 0) then
                read(u, '(A)') line ! Skip LOOKUP_TABLE line
                do i = 1, n; read(u, *) ids(i); end do
                print *, "DEBUG: Successfully read volume CellEntityIDS."
                return
            end if
        end do
400     print *, "WARNING: CellEntityIDS not found in source, defaulting to 1."
    end subroutine

    function label_bc(mx, my, x0, x1, y0, y1) result(id)
        real(dp), intent(in) :: mx, my, x0, x1, y0, y1
        integer :: id
        real(dp) :: tol; tol = (x1 - x0) * 1e-3
        if      (abs(my - y0) < tol) then; id = 101 ! Bottom
        else if (abs(mx - x1) < tol) then; id = 101 ! Right
        else if (abs(my - y1) < tol) then; id = 102 ! Top
        else if (abs(mx - x0) < tol) then; id = 102 ! Left
        else; id = 0; endif
    end function

    subroutine parse_points(u, p, n)
        integer, intent(in) :: u
        real(dp), allocatable, intent(out) :: p(:,:)
        integer, intent(out) :: n; integer :: i; character(len=256) :: line
        rewind(u)
        do while (.true.)
            read(u, '(A)', end=100) line
            if (index(line, "POINTS") > 0) exit
        end do
        read(line(index(line, "POINTS")+7:), *) n
        allocate(p(n, 2))
        do i = 1, n; read(u, *) p(i, 1), p(i, 2); end do
        return
100     stop "ERROR: POINTS not found"
    end subroutine

    subroutine parse_cells(u, corners, n, sz)
        integer, intent(in) :: u
        integer, allocatable, intent(out) :: corners(:,:)
        integer, intent(out) :: n, sz; integer :: i, n_node, ios; character(len=256) :: line
        do while (.true.)
            read(u, '(A)', end=200) line
            if (index(line, "CELLS") > 0 .and. index(line, "CELLS_") == 0) exit
        end do
        read(line(index(line, "CELLS")+6:), *) n, sz
        allocate(corners(n, 4))
        do i = 1, n
            read(u, *, iostat=ios) n_node, corners(i, 1), corners(i, 2), corners(i, 3), corners(i, 4)
            corners(i, :) = corners(i, :) + 1 
        end do
        return
200     stop "ERROR: CELLS not found"
    end subroutine

    function sort_pair(a, b) result(res)
        integer, intent(in) :: a, b; integer :: res(2)
        if (a < b) then; res = [a, b]; else; res = [b, a]; end if
    end function

    recursive subroutine quicksort_edges(a, first, last)
        type(t_edge), intent(inout) :: a(:); integer, intent(in) :: first, last
        integer :: i, j; type(t_edge) :: pivot, temp
        if (first >= last) return
        i = first; j = last; pivot = a((first + last) / 2)
        do while (i <= j)
            do while (a(i)%n(1) < pivot%n(1) .or. (a(i)%n(1) == pivot%n(1) .and. a(i)%n(2) < pivot%n(2))); i = i + 1; end do
            do while (a(j)%n(1) > pivot%n(1) .or. (a(j)%n(1) == pivot%n(1) .and. a(j)%n(2) > pivot%n(2))); j = j - 1; end do
            if (i <= j) then; temp = a(i); a(i) = a(j); a(j) = temp; i = i + 1; j = j - 1; end if
        end do
        if (first < j) call quicksort_edges(a, first, j)
        if (i < last) call quicksort_edges(a, i, last)
    end subroutine

subroutine write_vtk(f, p, q, m, b)
        character(len=*), intent(in) :: f
        real(dp), intent(in) :: p(:,:)
        integer, intent(in) :: q(:,:)
        integer, intent(in) :: m(:)
        type(t_edge), intent(in) :: b(:)
        integer :: i, u, nq, nb, tc
        
        nq = size(q, 1)
        nb = size(b, 1)
        tc = nq + nb
        
        open(newunit=u, file=trim(f), status='replace', action='write')
        
        write(u,'(A)') "# vtk DataFile Version 2.0"
        write(u,'(A)') "pincell_v1.0, Processed"
        write(u,'(A)') "ASCII"
        write(u,'(A)') "DATASET UNSTRUCTURED_GRID"
        
        write(u,'(A,I0,A)') "POINTS ", size(p, 1), " double"
        do i = 1, size(p, 1)
            write(u,*) p(i,1), p(i,2), 0.0_dp
        end do
        
        write(u,*) 
        write(u,'(A,I0,1X,I0)') "CELLS ", tc, (nq*5 + nb*3)
        do i = 1, nq
            write(u,'(I1,4(1X,I0))') 4, q(i,1)-1, q(i,2)-1, q(i,3)-1, q(i,4)-1
        end do
        do i = 1, nb
            write(u,'(I1,2(1X,I0))') 2, b(i)%n(1)-1, b(i)%n(2)-1
        end do
        
        write(u,*) 
        write(u,'(A,I0)') "CELL_TYPES ", tc
        do i = 1, nq
            write(u,'(I1)') 9 
        end do
        do i = 1, nb
            write(u,'(I1)') 3 
        end do
        
        write(u,*) 
        write(u,'(A,I0)') "CELL_DATA ", tc
        write(u,'(A)') "SCALARS CellEntityIds int 1"
        write(u,'(A)') "LOOKUP_TABLE default"
        
        do i = 1, nq
            write(u,'(I0)') m(i)
        end do
        do i = 1, nb
            write(u,'(I0)') b(i)%id
        end do
        
        close(u)
    end subroutine
end program VTK_NthOrder_Boundary