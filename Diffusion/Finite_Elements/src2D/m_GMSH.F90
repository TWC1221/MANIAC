module m_GMSH
    use m_constants
    use m_types, only: t_mesh
    implicit none
    private
    public :: parse

contains

    subroutine parse(filename, mesh, write_files)
        character(len=*), intent(in)  :: filename
        type(t_mesh),   intent(out)   :: mesh
        logical,          intent(in), optional :: write_files

        integer :: unit, ios
        integer :: nnode, nelem, total_ints
        integer :: i, j, k, npe, start, max_npe
        character(len=256) :: line, key, dtype
        logical :: do_write

        real(dp), allocatable :: nodes3(:,:)        
        integer,  allocatable :: connectivity(:)    
        integer,  allocatable :: elem_ptr(:)        
        integer,  allocatable :: elem_type(:)       
        integer,  allocatable :: material_all(:)    

        integer :: n_quads, n_lines
        integer :: iq, il

        do_write = .false.; if (present(write_files)) do_write = write_files
        
        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) stop "m_GMSH: cannot open file"

        ! --- POINTS section ---
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) stop "m_GMSH: POINTS section not found"
            if (index(line, "POINTS") > 0) exit
        end do
        read(line,*) key, nnode, dtype
        allocate(nodes3(3,nnode))
        do i = 1, nnode
            read(unit,*) nodes3(:,i)
        end do

        ! --- CELLS section ---
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) stop "m_GMSH: CELLS section not found"
            if (index(line, "CELLS") > 0) exit
        end do
        read(line,*) key, nelem, total_ints
        allocate(connectivity(total_ints))
        allocate(elem_ptr(nelem+1))
        read(unit,*) connectivity

        k = 1
        elem_ptr(1) = 1
        max_npe = 0
        do i = 1, nelem
            npe = connectivity(k)
            max_npe = max(max_npe, npe)
            elem_ptr(i+1) = elem_ptr(i) + 1 + npe
            k = elem_ptr(i+1)
        end do

        ! --- CELL_TYPES section ---
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) stop "m_GMSH: CELL_TYPES section not found"
            if (index(line, "CELL_TYPES") > 0) exit
        end do
        allocate(elem_type(nelem))
        do i = 1, nelem
            read(unit,*) elem_type(i)
        end do

        ! --- CELL_DATA (Materials) ---
        allocate(material_all(nelem)); material_all = 1
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit 
            if (index(line, "CELL_DATA") > 0) exit
        end do

        if (ios == 0) then
            read(unit,*) ! SCALARS...
            read(unit,*) ! LOOKUP_TABLE...
            do i = 1, nelem
                read(unit,*,iostat=ios) material_all(i)
                if (ios /= 0) material_all(i) = 1
            end do
        end if
        close(unit)

        mesh%dim     = 2
        mesh%n_nodes = nnode
        
        n_quads = count(elem_type == 9 .or. elem_type == 28)
        n_lines = count(elem_type == 3 .or. elem_type == 21)
        
        mesh%nloc = 0
        if (any(elem_type == 28)) then
            mesh%nloc = 9
        else if (any(elem_type == 9)) then
            mesh%nloc = 4
        end if

        allocate(mesh%nodes(mesh%n_nodes, 2))
        mesh%nodes(:,1) = nodes3(1,:)   
        mesh%nodes(:,2) = nodes3(2,:)   

        allocate(mesh%elems(n_quads, mesh%nloc))
        allocate(mesh%edges(n_lines, 3))
        allocate(mesh%mats(n_quads))
        allocate(mesh%edge_mats(n_lines))  

        mesh%edges = 0
        mesh%edge_mats = 0     

        iq = 0; il = 0
        do i = 1, nelem
            start = elem_ptr(i)
            npe   = connectivity(start)

            select case (elem_type(i))
            case (9, 28)
                iq = iq + 1
                mesh%mats(iq) = material_all(i)
                do j = 1, npe
                    mesh%elems(iq, j) = connectivity(start + j) + 1
                end do
            case (3, 21)
                il = il + 1
                mesh%edge_mats(il) = material_all(i) 
                do j = 1, npe
                    mesh%edges(il, j) = connectivity(start + j) + 1
                end do
            end select
        end do

        mesh%n_elems = n_quads
        mesh%n_edges = n_lines

        if (do_write) then
            call write_nodes(nodes3)
            call write_elements(nelem, connectivity, elem_ptr, elem_type)
            call write_materials(material_all)
            call write_edges(n_lines, mesh%edges, mesh%edge_mats)
            print *, "Mesh successfully parsed and written to .dat files."
        end if

        if (allocated(nodes3))        deallocate(nodes3)
        if (allocated(connectivity))  deallocate(connectivity)
        if (allocated(elem_ptr))      deallocate(elem_ptr)
        if (allocated(elem_type))     deallocate(elem_type)
        if (allocated(material_all))  deallocate(material_all)
    end subroutine parse

    subroutine write_edges(n_lines, edges, edge_mats)
        integer, intent(in) :: n_lines, edges(:,:), edge_mats(:)
        integer :: i, jj, unit
        open(newunit=unit, file="edges.dat", status="replace")
        write(unit,'(A)') '# id     material        nodes...'
        do i = 1, n_lines
            write(unit,'(I8,I12)',advance="no") i, edge_mats(i)
            do jj = 1, size(edges, 2)
                if (edges(i,jj) > 0) write(unit,'(I8)',advance="no") edges(i,jj)
            end do
            write(unit,*)
        end do
        close(unit)
    end subroutine write_edges

    subroutine write_nodes(nodes3)
        real(dp), intent(in) :: nodes3(:,:)
        integer :: i, unit
        open(newunit=unit, file="nodes.dat", status="replace")
        write(unit,'(A)') '# id          x                      y                      z'
        do i = 1, size(nodes3,2)
            write(unit,'(I8,3ES25.15)') i, nodes3(:,i)
        end do
        close(unit)
    end subroutine write_nodes

    subroutine write_elements(nelem, connectivity, elem_ptr, elem_type)
        integer, intent(in) :: nelem, connectivity(:), elem_ptr(:), elem_type(:)
        integer :: i, jj, npe, start, unit
        open(newunit=unit, file="elements.dat", status="replace")
        write(unit,'(A)') '# id     vtk_type        nodes...'
        do i = 1, nelem
            start = elem_ptr(i)
            npe   = connectivity(start)
            write(unit,'(I8,I12)',advance="no") i, elem_type(i)
            do jj = 1, npe
                write(unit,'(I8)',advance="no") connectivity(start+jj) + 1
            end do
            write(unit,*)
        end do
        close(unit)
    end subroutine write_elements

    subroutine write_materials(material_all)
        integer, intent(in) :: material_all(:)
        integer :: i, unit
        open(newunit=unit, file="mats.dat", status="replace")
        write(unit,'(A)') '# element_id   material'
        do i = 1, size(material_all)
            write(unit,'(I8, I12)') i, material_all(i)
        end do
        close(unit)
    end subroutine write_materials

end module m_GMSH