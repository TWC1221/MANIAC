module m_material
  use m_constants
  use m_types
  implicit none
  contains

  subroutine InitialiseMaterials(mats, mesh, n_groups, filename, printout)
      type(t_material), allocatable, intent(out) :: mats(:)
      type(t_mesh), intent(in) :: mesh
      integer, intent(in) :: n_groups
      character(len=*), intent(in) :: filename
      logical, intent(in) :: printout
      integer, allocatable :: unique_ids(:)
      integer :: i, max_id, n_unique

      max_id = max(maxval(mesh%mats), maxval(mesh%edge_mats))
      allocate(mats(max_id))

      allocate(unique_ids(max_id))
      n_unique = 0
      do i = 1, max_id
          if (any(mesh%mats == i) .or. any(mesh%edge_mats == i)) then
              n_unique = n_unique + 1
              unique_ids(n_unique) = i
          end if
      end do

      call ParseMaterialDeck(mats, unique_ids(1:n_unique), n_groups, filename)

      do i = 1, max_id
          if (allocated(mats(i)%SigA)) call UpdateComputables(mats(i), n_groups)
      end do

      if (printout) call PrintMaterialSummary(mats, n_groups)
  end subroutine InitialiseMaterials

  subroutine ParseMaterialDeck(mats, requested_ids, n_groups, filename)
      type(t_material), intent(inout) :: mats(:)
      integer, intent(in)             :: requested_ids(:), n_groups
      character(len=*), intent(in)    :: filename
      integer :: u, id, g, ios, i
      character(len=1024) :: line
      
      open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      
      do
          read(u, '(A)', iostat=ios) line
          if (ios /= 0) exit
          line = adjustl(line)
          
          if (line == '' .or. line(1:1) == '-' .or. index(line, 'matID') > 0) cycle

          read(line, *, iostat=ios) id
          if (ios /= 0) cycle
          
          if (.not. any(requested_ids == id)) cycle

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
              
              if (g == 1) then
                  read(line, *, iostat=ios) id, mats(id)%SigA(g), mats(id)%Nu(g), &
                                mats(id)%SigF(g), mats(id)%Chi(g),    &
                                mats(id)%Src(g), (mats(id)%SigmaS(g, i), i=1, n_groups), &
                                mats(id)%name
              else
                  read(line, *, iostat=ios) id, mats(id)%SigA(g), mats(id)%Nu(g), &
                                mats(id)%SigF(g), mats(id)%Chi(g),    &
                                mats(id)%Src(g), (mats(id)%SigmaS(g, i), i=1, n_groups)
              end if
          end do
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

subroutine PrintMaterialSummary(mats, n_groups)
      type(t_material), intent(in) :: mats(:)
      integer, intent(in) :: n_groups
      integer :: i, g, gp
      character(len=256) :: fmt_data
      character(len=170) :: line, sep
      integer :: total_width

      total_width = 61 + (n_groups * 13) 
      
      line = " |" // repeat('=', total_width-3) // "|"
      sep  = " |-----+------------+------------+------------+------------+" // repeat('-', total_width-61) // "|"

      write(fmt_data, '(A, I2, A, A)') &
        '(A, I3,X, 4(" | ", ES10.3), " | ", ', n_groups, '(ES11.4, " |"), A)'

      write(*,'(A)') line
      write(*,'(A)') " |" // repeat(' ', (total_width-28)/2) // "MATERIAL DATABASE SUMMARY" // &
                     repeat(' ', (total_width-28)/2) // "|"
      write(*,'(A)') line

      do i = 1, size(mats)
          if (.not. allocated(mats(i)%SigA)) cycle

          write(*,'(A, I3, A, A, A, T' // to_str(total_width) // ', A)') &
            " | [M] Material ID: ", i, " (", trim(mats(i)%name), ")", "|"
          write(*,'(A)') " |" // repeat('-', total_width-3) // "|"
          
          write(*,'(A, T63, A, T' // to_str(total_width) // ', A)') &
            " | Grp |    SigT    |    SigA    |   NuSigF   |     Chi    |", &
            "Scattering Matrix (From Row \ To Column)", "|"
          write(*,'(A)') sep

          do g = 1, n_groups
              write(*, fmt_data) &
                  " |", g, mats(i)%SigmaT(g), mats(i)%SigA(g), mats(i)%NuSigF(g), mats(i)%Chi(g), &
                  (mats(i)%SigmaS(g, gp), gp=1, n_groups), ""
          end do

          write(*,'(A)') line
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