module m_GMSH
    use omp_lib
    use m_constants
    use m_types, only: t_mesh
    use m_mesh_quality, only: check_mesh_quality
    implicit none
    private
    public :: parse_GMSH

contains
      
    pure function derive_case_nametag(filename) result(tag)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: tag
        integer :: pos

        pos = index(filename, '/', back=.true.)
        if (pos == 0) pos = index(filename, '\', back=.true.)

        tag = filename(pos+1:)
    end function derive_case_nametag

    subroutine parse_GMSH(filename, mesh, write_files, verbose)
        character(len=*), intent(in)  :: filename
        type(t_mesh),   intent(out)   :: mesh
        logical,          intent(in), optional :: write_files
        logical,          intent(in), optional :: verbose
        logical, save :: first_call = .true.

        integer :: unit, ios
        integer :: nnode, nelem, total_ints
        integer :: i, j, k, npe, start, max_npe
        character(len=256) :: line, key, dtype
        logical :: do_write
        logical :: print_log

        integer, allocatable :: temp_elements_nodes(:,:)
        real(dp), allocatable :: nodes3(:,:)        
        integer,  allocatable :: connectivity(:)    
        integer,  allocatable :: elem_ptr(:)        
        integer,  allocatable :: elem_type(:)       
        integer,  allocatable :: material_all(:)    

        integer :: n_quads, n_lines, n_edge_nodes
        integer :: iq, il

        do_write = .false.; if (present(write_files)) do_write = write_files
        print_log = .false.; if (present(verbose)) print_log = verbose
        
        if (print_log) write(*,'(A)') " [ SYSTEM ] :: Loading Mesh: " // derive_case_nametag(trim(filename))
        if (print_log) write(*,'(A, I0)') " [ SYSTEM ] :: OMP_NUM_THREADS: ", omp_get_max_threads()

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

        mesh%dim     = 3
        mesh%n_nodes = nnode
        
        n_quads = count(elem_type == 12) ! VTK_HEXAHEDRON
        n_lines = count(elem_type == 9)  ! VTK_QUAD
        
        mesh%nloc = 0
        n_edge_nodes = 0

        do i = 1, nelem
            if (elem_type(i) == 12) then
                mesh%nloc = max(mesh%nloc, connectivity(elem_ptr(i)))
            else if (elem_type(i) == 9) then
                n_edge_nodes = max(n_edge_nodes, connectivity(elem_ptr(i)))
            end if
        end do

        allocate(mesh%nodes(mesh%n_nodes, 3))
        mesh%nodes(:,1) = nodes3(1,:)   
        mesh%nodes(:,2) = nodes3(2,:)   
        mesh%nodes(:,3) = nodes3(3,:)

        allocate(mesh%elems(n_quads, mesh%nloc))
        allocate(mesh%edges(n_lines, n_edge_nodes))
        allocate(mesh%mats(n_quads))
        allocate(mesh%edge_mats(n_lines))  

        mesh%elems = 0
        mesh%edges = 0
        mesh%edge_mats = 0     

        iq = 0; il = 0
        do i = 1, nelem
            start = elem_ptr(i)
            npe   = connectivity(start)

            select case (elem_type(i))
            case (12)
                iq = iq + 1
                mesh%mats(iq) = material_all(i)
                do j = 1, npe
                    mesh%elems(iq, j) = connectivity(start + j) + 1
                end do
            case (9)
                il = il + 1
                mesh%edge_mats(il) = material_all(i) 
                do j = 1, npe
                    mesh%edges(il, j) = connectivity(start + j) + 1
                end do
            end select
        end do

        mesh%n_elems = n_quads
        mesh%n_edges = n_lines

        if (first_call) then
            allocate(temp_elements_nodes(max_npe, nelem))
            temp_elements_nodes = 0
            do i = 1, nelem
                start = elem_ptr(i)
                npe = connectivity(start)
                do j = 1, npe
                    temp_elements_nodes(j, i) = connectivity(start + j) + 1
                end do
            end do

            call check_mesh_quality(nodes3, temp_elements_nodes, elem_type, nelem, bc_ids=material_all, printout=print_log)
            deallocate(temp_elements_nodes)
            first_call = .false.
        end if

        if (do_write) then
            call write_nodes(nodes3)
            call write_elements(nelem, connectivity, elem_ptr, elem_type)
            call write_edges(n_lines, mesh%edges, mesh%edge_mats)
            call write_materials(mesh%n_elems, mesh%mats)
            if (print_log) write(*,'(A)') " [  I/O   ] :: MESH PARSED: .dat files WRITTEN"
        end if

        if (allocated(nodes3))        deallocate(nodes3)
        if (allocated(connectivity))  deallocate(connectivity)
        if (allocated(elem_ptr))      deallocate(elem_ptr)
        if (allocated(elem_type))     deallocate(elem_type)
        if (allocated(material_all))  deallocate(material_all)
    end subroutine parse_GMSH

    subroutine write_edges(n_lines, edges, edge_mats)
        integer, intent(in) :: n_lines, edges(:,:), edge_mats(:)
        integer :: i, jj, unit
        open(newunit=unit, file="../output/edges.dat", status="replace")
        write(unit,'(A)') '# Edge_ID    BC_ID      Nodes...'
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
        open(newunit=unit, file="../output/nodes.dat", status="replace")
        write(unit,'(A)') '# id          x                      y                      z'
        do i = 1, size(nodes3,2)
            write(unit,'(I8,3ES25.15)') i, nodes3(:,i)
        end do
        close(unit)
    end subroutine write_nodes

    subroutine write_elements(nelem, connectivity, elem_ptr, elem_type)
        integer, intent(in) :: nelem, connectivity(:), elem_ptr(:), elem_type(:)
        integer :: i, jj, npe, start, unit
        open(newunit=unit, file="../output/elements.dat", status="replace")
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

    subroutine write_materials(n_elems, mats)
        integer, intent(in) :: n_elems, mats(:)
        integer :: i, unit
        open(newunit=unit, file="../output/mats.dat", status="replace")
        write(unit,'(A)') '# Elem_ID    Mat_ID'
        do i = 1, n_elems
            write(unit,'(I8, I10)') i, mats(i)
        end do
        close(unit)
    end subroutine write_materials

end module m_GMSH