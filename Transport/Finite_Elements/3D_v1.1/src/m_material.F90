module m_material
  use m_constants
  use m_types, only: t_material, t_mesh, timer_start, timer_stop, int_to_str
  implicit none
  private
  public :: InitialiseMaterials, ParseMaterialDeck, UpdateComputables, &
            PrintMaterialSummary

  contains

  subroutine InitialiseMaterials(mats, mesh, n_groups, filename, printout, file_unit)
      type(t_material), allocatable, intent(out) :: mats(:)
      type(t_mesh), intent(in) :: mesh
      integer, intent(in) :: n_groups
      character(len=*), intent(in) :: filename
      integer, intent(in) :: file_unit ! Always writes to this unit
      logical, intent(in) :: printout
      integer, allocatable :: unique_ids(:)
      integer :: i, max_id, n_unique

      call timer_start('I/O  - Parse Materials')

      max_id = max(maxval(mesh%material_ids), maxval(mesh%boundary_ids))
      allocate(mats(max_id))

      allocate(unique_ids(max_id))
      n_unique = 0
      do i = 1, max_id
          if (any(mesh%material_ids == i) .or. any(mesh%boundary_ids == i)) then
              n_unique = n_unique + 1
              unique_ids(n_unique) = i
          end if
      end do

      call ParseMaterialDeck(mats, unique_ids(1:n_unique), n_groups, filename)

      do i = 1, max_id
          if (allocated(mats(i)%SigA)) call UpdateComputables(mats(i), n_groups)
      end do

      call PrintMaterialSummary(mats, n_groups, file_unit, printout)

      call timer_stop('I/O - Parse Materials')
  end subroutine InitialiseMaterials

  subroutine ParseMaterialDeck(mats, requested_ids, n_groups, filename)
      type(t_material), intent(inout) :: mats(:)
      integer, intent(in)             :: requested_ids(:), n_groups
      character(len=*), intent(in)    :: filename
      integer :: u, id, g, ios, i, max_id
      character(len=1024) :: line
      integer, parameter :: N_GROUPS_IN_FILE = 7
      real(dp) :: temp_scat(N_GROUPS_IN_FILE)
      logical, allocatable :: processed_ids(:)
      
      open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      max_id = size(mats)
      allocate(processed_ids(max_id)); processed_ids = .false.
      
      do
          read(u, '(A)', iostat=ios) line
          if (ios /= 0) exit
          line = adjustl(line)
          
          if (line == '' .or. line(1:1) == '-' .or. index(line, 'matID') > 0) cycle

          read(line, *, iostat=ios) id
          if (ios /= 0) cycle
          
          if (id > max_id .or. id < 1) cycle
          if (.not. any(requested_ids == id) .or. processed_ids(id)) cycle

          if (.not. allocated(mats(id)%SigA)) then
              allocate(mats(id)%SigA(n_groups), mats(id)%Nu(n_groups),    &
                       mats(id)%SigF(n_groups), mats(id)%NuSigF(n_groups),&
                       mats(id)%Chi(n_groups),  mats(id)%Src(n_groups),   &
                       mats(id)%SigmaS(n_groups, n_groups),               &
                       mats(id)%D(n_groups), mats(id)%SigmaT(n_groups),   &
                       mats(id)%SigmaR(n_groups))
              mats(id)%SigmaS = 0.0_dp
              mats(id)%name   = 'Unknown'
          end if

          do g = 1, n_groups
              if (g > 1) then
                  read(u, '(A)', iostat=ios) line
                  if (ios /= 0) then
                      write(*,*) "Error reading material file: unexpected end or failure at ID ", id, " Group ", g
                      exit
                  end if
              end if
              
              temp_scat = 0.0_dp
              if (g == 1) then
                  read(line, *, iostat=ios) id, mats(id)%SigA(g), mats(id)%Nu(g), &
                                mats(id)%SigF(g), mats(id)%Chi(g),    &
                                mats(id)%Src(g), (temp_scat(i), i=1, N_GROUPS_IN_FILE), &
                                mats(id)%name
              else
                  read(line, *, iostat=ios) id, mats(id)%SigA(g), mats(id)%Nu(g), &
                                mats(id)%SigF(g), mats(id)%Chi(g),    &
                                mats(id)%Src(g), (temp_scat(i), i=1, N_GROUPS_IN_FILE)
              end if

              if (ios /= 0) then
                  write(*,*) "Error parsing line for material ", id, " group ", g
                  write(*,*) "Line: ", trim(line)
                  stop "Material parsing failed."
              end if
              mats(id)%SigmaS(g, 1:n_groups) = temp_scat(1:n_groups)
          end do
          processed_ids(id) = .true.
      end do
      close(u)
  end subroutine ParseMaterialDeck

  subroutine UpdateComputables(mat, n_groups)
      type(t_material), intent(inout) :: mat
      integer, intent(in) :: n_groups
      integer :: g

      if (.not. allocated(mat%SigmaT)) allocate(mat%SigmaT(n_groups))
      if (.not. allocated(mat%SigmaR)) allocate(mat%SigmaR(n_groups))
      if (.not. allocated(mat%D))      allocate(mat%D(n_groups))

      if (sum(mat%Chi) > 1e-12_dp) mat%Chi = mat%Chi / sum(mat%Chi)

      do g = 1, n_groups
          mat%NuSigF(g) = mat%Nu(g) * mat%SigF(g)
          mat%SigmaT(g) = mat%SigA(g) + sum(mat%SigmaS(g, :))
          mat%SigmaR(g) = mat%SigmaT(g) - mat%SigmaS(g, g)
          
          if (mat%SigmaT(g) > 1e-10_dp) then
              mat%D(g) = 1.0_dp / (3.0_dp * mat%SigmaT(g))
          else
              print*,"ERROR: D"
          end if
      end do
  end subroutine UpdateComputables

subroutine PrintMaterialSummary(mats, n_groups, file_unit, printout)
      type(t_material), intent(in) :: mats(:)
      integer, intent(in) :: n_groups
      integer, intent(in) :: file_unit ! Always writes to this unit
      integer :: i, g, gp, u, n_out, u_list(2)
      logical, intent(in) :: printout 
      character(len=256) :: fmt_data
      character(len=170) :: line, sep
      integer :: total_width

      total_width = 61 + (n_groups * 13) 
      
      line = " |" // repeat('=', total_width-3) // "|"
      sep  = " |-----+------------+------------+------------+------------+" // repeat('-', total_width-61) // "|"
      
      u_list(1) = file_unit
      n_out = 1
      if (printout) then
          u_list(2) = 6
          n_out = 2
      end if

      write(fmt_data, '(A, I2, A, A)') &
        '(A, I3,X, 4(" | ", ES10.3), " | ", ', n_groups, '(ES11.4, " |"), A)'

      do u = 1, n_out
          write(u_list(u),'(A)') line
          write(u_list(u),'(A)') " |" // repeat(' ', (total_width-28)/2) // "MATERIAL DATABASE SUMMARY" // &
                         repeat(' ', (total_width-28)/2) // "|"
          write(u_list(u),'(A)') line
      end do

      do i = 1, size(mats)
          if (.not. allocated(mats(i)%SigA)) cycle

          do u = 1, n_out
              write(u_list(u),'(A, I3, A, A, A, T' // to_str(total_width) // ', A)') &
                " | [M] Material ID: ", i, " (", trim(mats(i)%name), ")", "|"
              write(u_list(u),'(A)') " |" // repeat('-', total_width-3) // "|"
              
              write(u_list(u),'(A, T63, A, T' // to_str(total_width) // ', A)') &
                " | Grp |    SigT    |    SigA    |   NuSigF   |     Chi    |", &
                "Scattering Matrix (From Row \ To Column)", "|"
              write(u_list(u),'(A)') sep

              do g = 1, n_groups
                  write(u_list(u), fmt_data) &
                      " |", g, mats(i)%SigmaT(g), mats(i)%SigA(g), mats(i)%NuSigF(g), mats(i)%Chi(g), &
                      (mats(i)%SigmaS(g, gp), gp=1, n_groups), ""
              end do
              write(u_list(u),'(A)') line
          end do
      end do

  contains
      function to_str(n) result(str)
          integer, intent(in) :: n
          character(len=4) :: str
          write(str, '(I4)') n
          str = adjustl(str)
      end function to_str
  end subroutine PrintMaterialSummary

end module m_material