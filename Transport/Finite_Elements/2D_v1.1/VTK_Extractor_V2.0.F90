program VTK_Nth_Extractor
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    
    type t_edge
        integer :: n(10)      
        integer :: n_nodes    
        integer :: sort_key(2)
        integer :: id         
    end type t_edge

    real(dp), allocatable :: pts(:,:)
    integer, allocatable  :: all_nodes(:,:), cell_types(:)
    integer, allocatable  :: cell_entity_ids(:), pin_ids(:)
    type(t_edge), allocatable :: all_e(:), b_e(:)
    
    real(dp) :: xmin, xmax, ymin, ymax, mx, my
    integer  :: npts, ncells, n_b, i, j, k, npc, total_sz, order, npe
    character(len=100) :: infile, outfile, filename
    logical :: exists

    filename = "C5G7vol_2_5_f.vtk"
    infile = "../input/raw/" // trim(filename)
    outfile = "../input/" // trim(filename)
    print*,"Input File:" // filename

    inquire(file=trim(infile), exist=exists)
    if (.not. exists) stop "ERROR: Input file not found"

    open(10, file=trim(infile), status='old', action='read')
    call parse_points(10, pts, npts)
    call parse_cells_flexible(10, all_nodes, ncells, npc, total_sz)
    call parse_cell_types(10, cell_types, ncells)
    call parse_scalars(10, cell_entity_ids, ncells, "CellEntityIds") 
    call parse_scalars(10, pin_ids, ncells, "PinID") 
    close(10)

    order = nint(sqrt(real(npc, dp))) - 1
    npe   = order + 1
    
    xmin = minval(pts(:,1)); xmax = maxval(pts(:,1))
    ymin = minval(pts(:,2)); ymax = maxval(pts(:,2))
    
    allocate(all_e(ncells * 4 + 10))
    k = 0
    do i = 1, ncells
        do j = 1, 4
            k = k + 1
            call extract_vtk_lagrange_face(all_nodes(i,:), j, order, all_e(k))
        end do
    end do

    call quicksort_edges(all_e, 1, k)
    allocate(b_e(k))
    n_b = 0; i = 1
    do while (i <= k)
        if (i < k .and. all(all_e(i)%sort_key == all_e(i+1)%sort_key)) then
            i = i + 2 
        else
            n_b = n_b + 1
            b_e(n_b) = all_e(i)
            mx = (pts(b_e(n_b)%n(1), 1) + pts(b_e(n_b)%n(2), 1)) / 2.0_dp
            my = (pts(b_e(n_b)%n(1), 2) + pts(b_e(n_b)%n(2), 2)) / 2.0_dp
            b_e(n_b)%id = label_bc(mx, my, xmin, xmax, ymin, ymax)
            i = i + 1
        end if
    end do

    call write_vtk_final(outfile, pts, all_nodes, cell_types, &
                         cell_entity_ids, pin_ids, b_e(1:n_b), npc, npe, order)
    
contains

    subroutine extract_vtk_lagrange_face(nodes, side, p, edge)
        integer, intent(in) :: nodes(:), side, p
        type(t_edge), intent(out) :: edge
        integer :: npe, m_idx
        npe = p + 1; edge%n_nodes = npe; m_idx = 5
        select case(side)
        case(1) ! Bottom
            edge%n(1)=nodes(1); edge%n(2)=nodes(2)
            if(p>1) edge%n(3:npe)=nodes(m_idx : m_idx+(p-2))
        case(2) ! Right
            edge%n(1)=nodes(2); edge%n(2)=nodes(3)
            if(p>1) edge%n(3:npe)=nodes(m_idx+(p-1) : m_idx+2*(p-1)-1)
        case(3) ! Top
            edge%n(1)=nodes(3); edge%n(2)=nodes(4)
            if(p>1) edge%n(3:npe)=nodes(m_idx+2*(p-1) : m_idx+3*(p-1)-1)
        case(4) ! Left
            edge%n(1)=nodes(4); edge%n(2)=nodes(1)
            if(p>1) edge%n(3:npe)=nodes(m_idx+3*(p-1) : m_idx+4*(p-1)-1)
        end select
        edge%sort_key(1)=min(edge%n(1), edge%n(2))
        edge%sort_key(2)=max(edge%n(1), edge%n(2))
    end subroutine

    subroutine write_vtk_final(f, pts, nodes, types, m1, m2, b, npc, npe, order)
        character(len=*), intent(in) :: f; real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: nodes(:,:), types(:), m1(:), m2(:), npc, npe, order
        type(t_edge), intent(in) :: b(:)
        integer :: i, u, nq, nb, tc, sz_calc, vtk_line
        nq=size(nodes,1); nb=size(b,1); tc=nq+nb
        sz_calc=(nq*(npc+1))+(nb*(npe+1))
        vtk_line = 3; if(order==2) vtk_line=21; if(order==3) vtk_line=35
        open(newunit=u, file=trim(f), status='replace')
        write(u,'(A)') "# vtk DataFile Version 2.0"; write(u,'(A)') "NTH Order Polynomial Mesh"
        write(u,'(A)') "ASCII"; write(u,'(A)') "DATASET UNSTRUCTURED_GRID"
        write(u,'(A,I0,A)') "POINTS ", size(pts,1), " double"
        do i=1,size(pts,1); write(u,'(3(ES16.8,1X))') pts(i,1),pts(i,2),0.0_dp; end do
        write(u,*) ""
        write(u,'(A,I0,1X,I0)') "CELLS ", tc, sz_calc
        do i=1,nq; write(u,'(I0,99(1X,I0))') npc, nodes(i,:)-1; end do
        do i=1,nb; write(u,'(I0,99(1X,I0))') npe, b(i)%n(1:npe)-1; end do
        write(u,*) ""
        write(u,'(A,I0)') "CELL_TYPES ", tc
        do i=1,nq; write(u,'(I0)') types(i); end do
        do i=1,nb; write(u,'(I0)') vtk_line; end do
        write(u,*) ""
        write(u,'(A,I0)') "CELL_DATA ", tc
        write(u,'(A)') "SCALARS CellEntityIds int 1"
        write(u,'(A)') "LOOKUP_TABLE default"
        do i=1,nq; write(u,'(I0)') m1(i); end do
        do i=1,nb; write(u,'(I0)') b(i)%id; end do
        write(u,*) ""
        write(u,'(A)') "SCALARS PinID int 1"
        write(u,'(A)') "LOOKUP_TABLE default"
        do i=1,nq; write(u,'(I0)') m2(i); end do
        do i=1,nb; write(u,'(I0)') 0; end do
        close(u)
    end subroutine

    subroutine parse_points(u, p, n)
        integer, intent(in) :: u; integer, intent(out) :: n; real(dp), allocatable, intent(out) :: p(:,:)
        character(len=256) :: ln; integer :: i; rewind(u)
        do; read(u,'(A)') ln; if(index(ln,"POINTS")>0) exit; end do
        read(ln(7:),*) n; allocate(p(n,2))
        do i=1,n; read(u,*) p(i,1), p(i,2); end do
    end subroutine

    subroutine parse_cells_flexible(u, nds, n, npc, sz)
        integer, intent(in) :: u; integer, allocatable, intent(out) :: nds(:,:)
        integer, intent(out) :: n, npc, sz; character(len=256) :: ln; integer :: i, d
        rewind(u); do; read(u,'(A)') ln; if(index(ln,"CELLS")>0.and.index(ln,"CELLS_")==0) exit; end do
        read(ln(6:),*) n, sz; read(u,*) npc; backspace(u); allocate(nds(n, npc))
        do i=1,n; read(u,*) d, nds(i,1:d); nds(i,:)=nds(i,:)+1; end do
    end subroutine

    subroutine parse_cell_types(u, t, n)
        integer, intent(in) :: u, n; integer, allocatable, intent(out) :: t(:)
        character(len=256) :: ln; integer :: i; rewind(u)
        do; read(u,'(A)') ln; if(index(ln,"CELL_TYPES")>0) exit; end do
        allocate(t(n)); do i=1,n; read(u,*) t(i); end do
    end subroutine

    subroutine parse_scalars(u, s, n, nm)
        integer, intent(in) :: u, n; character(len=*) :: nm; integer, allocatable :: s(:)
        character(len=256) :: ln; integer :: i; rewind(u); allocate(s(n)); s=0
        do; read(u,'(A)',end=10) ln; if(index(ln,"SCALARS")>0.and.index(ln,trim(nm))>0) then
        read(u,'(A)') ln; do i=1,n; read(u,*) s(i); end do; return; end if; end do
        10 continue
    end subroutine

    function label_bc(mx, my, x0, x1, y0, y1) result(id)
        real(dp), intent(in) :: mx, my, x0, x1, y0, y1; integer :: id; real(dp) :: tol
        tol = (x1 - x0) * 1e-3; id = 0
        if(abs(my-y0)<tol) id=101; if(abs(mx-x1)<tol) id=102
        if(abs(my-y1)<tol) id=103; if(abs(mx-x0)<tol) id=104
    end function

    recursive subroutine quicksort_edges(a, first, last)
        type(t_edge), intent(inout) :: a(:); integer, intent(in) :: first, last
        integer :: i, j; type(t_edge) :: pvt, tmp
        if (first >= last) return
        i=first; j=last; pvt=a((first+last)/2)
        do while (i <= j)
            do while (compare(a(i), pvt) < 0)
                i = i + 1
            end do
            do while (compare(a(j), pvt) > 0)
                j = j - 1
            end do
            if (i <= j) then
                tmp = a(i); a(i) = a(j); a(j) = tmp
                i = i + 1; j = j - 1
            end if
        end do
        if (first < j) call quicksort_edges(a, first, j)
        if (i < last) call quicksort_edges(a, i, last)
    end subroutine

    function compare(e1, e2) result(res)
        type(t_edge), intent(in) :: e1, e2; integer :: res
        if (e1%sort_key(1) < e2%sort_key(1)) then
            res = -1
        else if (e1%sort_key(1) > e2%sort_key(1)) then
            res = 1
        else if (e1%sort_key(2) < e2%sort_key(2)) then
            res = -1
        else if (e1%sort_key(2) > e2%sort_key(2)) then
            res = 1
        else
            res = 0
        end if
    end function
end program