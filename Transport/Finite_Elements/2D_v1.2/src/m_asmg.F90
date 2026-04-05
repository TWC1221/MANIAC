module m_asmg
    use m_constants
    use m_outpVTK_dfem
    use m_types
    implicit none
    private
    public :: read_asmg_mesh

contains

    subroutine read_asmg_mesh(filepath, mesh)
        character(len=*), intent(in) :: filepath
        type(t_mesh), intent(inout) :: mesh

        integer :: unit, iostatus, i, j, npe, pos
        character(len=1024) :: line
        integer :: v(4), v2(2)
        real(dp), allocatable :: temp_nodes3(:,:)
        integer :: p_count, patch_count, edge_count
        
        ! --- Pass 1: Count Elements and Nodes ---
        open(newunit=unit, file=filepath, status='old', action='read')
        patch_count = 0
        edge_count = 0
        p_count = 0
        do
            read(unit, '(A)', iostat=iostatus) line
            if (iostatus /= 0) exit
            line = adjustl(line)
            if (index(line, '$2D_Patch_Description_Start') > 0) then
                patch_count = patch_count + 1
            else if (index(line, '$1D_Patch_Description_Start') > 0) then
                edge_count = edge_count + 1
            else if (index(line, 'POINTS') == 1) then
                read(line(index(line,' '):), *) p_count
            end if
        end do

        mesh%n_elems = patch_count
        mesh%n_nodes = p_count
        mesh%n_edges = edge_count
        mesh%dim     = 2
        mesh%n_faces_per_elem = 4
        npe = 4 ! Quads

        allocate(mesh%nodes(mesh%n_nodes, 2))
        allocate(mesh%elems(mesh%n_elems, npe))
        allocate(mesh%material_ids(mesh%n_elems))
        allocate(mesh%boundary_ids(mesh%n_edges))
        allocate(mesh%edges(mesh%n_edges, 2))
        allocate(mesh%pin_ids(mesh%n_elems)); mesh%pin_ids = 0
        
        rewind(unit)

        ! --- Pass 2: Parse Data ---
        do
            read(unit, '(A)', iostat=iostatus) line
            if (iostatus /= 0) exit
            line = adjustl(line)

            if (index(line, 'POINTS') == 1) then
                i = 1
                do while (i <= mesh%n_nodes)
                    read(unit, *) mesh%nodes(i, 1), mesh%nodes(i, 2)
                    i = i + 1
                end do

            else if (index(line, '$$$ 1D_PATCHES_START $$$') > 0) then
                j = 1
                do while (j <= mesh%n_edges)
                    read(unit, '(A)') line
                    if (index(line, '$1D_Patch_Description_Start') > 0) then
                        do
                            read(unit, '(A)') line
                            line = adjustl(line)
                            if (index(line, 'BC:') == 1) then
                                call skip_to_val(unit, line)
                                read(line, *) mesh%boundary_ids(j)
                            else if (index(line, 'control_points:') == 1) then
                                call skip_to_val(unit, line)
                                read(line, *) v2
                                mesh%edges(j, 1) = v2(1) + 1
                                mesh%edges(j, 2) = v2(2) + 1
                            else if (index(line, '$1D_Patch_Description_End') > 0) then
                                exit
                            end if
                        end do
                        j = j + 1
                    end if
                end do

            else if (index(line, '$$$ 2D_PATCHES_START $$$') > 0) then
                j = 1
                do while (j <= mesh%n_elems)
                    read(unit, '(A)') line
                    if (index(line, '$2D_Patch_Description_Start') > 0) then
                        do
                            read(unit, '(A)') line
                            line = adjustl(line)
                            if (len_trim(line) == 0 .or. line(1:1) == '!') cycle

                            ! Robust same-line parsing for IDs
                            pos = index(line, 'Material_ID:')
                            if (pos > 0) then
                                read(line(pos+12:), *) mesh%material_ids(j)
                            end if

                            pos = index(line, 'Pin_ID:')
                            if (pos > 0) then
                                read(line(pos+7:), *) mesh%pin_ids(j)
                            end if

                            if (index(line, 'control_points:') == 1) then
                                call skip_to_val(unit, line)
                                read(line, *) v
                                mesh%elems(j, :) = v + 1
                            else if (index(line, '$2D_Patch_Description_End') > 0) then
                                exit
                            end if
                        end do
                        j = j + 1
                    end if
                end do
            end if
        end do
        close(unit)

        ! --- Writeouts for Debugging ---
        block
            character(len=512) :: out_dir
            out_dir = "../output/" // derive_case_nametag(filepath)
            call execute_command_line("mkdir -p " // trim(out_dir))
            
            allocate(temp_nodes3(mesh%n_nodes, 3))
            temp_nodes3(:,1:2) = mesh%nodes
            temp_nodes3(:,3)   = 0.0_dp

            call write_nodes(temp_nodes3, out_dir)
            call write_elements(mesh%n_elems, mesh%elems, out_dir)
            call write_edges(mesh%n_edges, mesh%edges, mesh%boundary_ids, out_dir)
            call write_materials(mesh%n_elems, mesh%material_ids, out_dir)
            call write_pins(mesh%n_elems, mesh%pin_ids, out_dir)
            
            deallocate(temp_nodes3)
        end block

        print *, "Successfully loaded ASMG mesh with boundaries: ", trim(filepath)
    end subroutine read_asmg_mesh

    subroutine skip_to_val(u, line)
        integer, intent(in) :: u
        character(len=*), intent(out) :: line
        do
            read(u, '(A)') line
            if (len_trim(line) > 0 .and. line(1:1) /= '!') exit
        end do
    end subroutine skip_to_val

    subroutine write_nodes(nodes3, out_dir)
        real(dp), intent(in) :: nodes3(:,:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/nodes.dat", status="replace")
        write(unit,'(A)') '# id          x                      y                      z'
        do i = 1, size(nodes3, 1)
            write(unit,'(I8,3ES25.15)') i, nodes3(i,:)
        end do
        close(unit)
    end subroutine write_nodes

    subroutine write_elements(nelem, elems, out_dir)
        integer, intent(in) :: nelem, elems(:,:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/elements.dat", status="replace")
        write(unit,'(A)') '# id     vtk_type        nodes...'
        do i = 1, nelem
            write(unit,'(I8,I12,4I8)') i, 9, elems(i, :) ! 9 = VTK_QUAD
        end do
        close(unit)
    end subroutine write_elements

    subroutine write_edges(n_edges, edges, bc_ids, out_dir)
        integer, intent(in) :: n_edges, edges(:,:), bc_ids(:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/edges.dat", status="replace")
        write(unit,'(A)') '# id     node1    node2    bc_id'
        do i = 1, n_edges
            write(unit,'(4I8)') i, edges(i,1), edges(i,2), bc_ids(i)
        end do
        close(unit)
    end subroutine write_edges

    subroutine write_materials(n_elems, mats, out_dir)
        integer, intent(in) :: n_elems, mats(:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/mats.dat", status="replace")
        write(unit,'(A)') '# Elem_ID    Mat_ID'
        do i = 1, n_elems
            write(unit,'(I8, I10)') i, mats(i)
        end do
        close(unit)
    end subroutine write_materials

    subroutine write_pins(n_elems, pins, out_dir)
        integer, intent(in) :: n_elems, pins(:)
        character(len=*), intent(in) :: out_dir
        integer :: i, unit
        open(newunit=unit, file=trim(out_dir)//"/pins.dat", status="replace")
        write(unit,'(A)') '# Elem_ID    Pin_ID'
        do i = 1, n_elems
            write(unit,'(I8, I10)') i, pins(i)
        end do
        close(unit)
    end subroutine write_pins
end module m_asmg
