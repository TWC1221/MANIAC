module m_materials
  use m_constants
  implicit none 

  type t_material
      character(len=32) :: name
      real(dp), allocatable :: D(:)       ! Diffusion coefficient [cm]
      real(dp), allocatable :: SigmaT(:)  ! Total Cross Section [cm^-1]
      real(dp), allocatable :: SigmaR(:)  ! Removal Cross Section [cm^-1]
      real(dp), allocatable :: SigA(:)    ! Absorption Cross Section [cm^-1]
      real(dp), allocatable :: Nu(:)
      real(dp), allocatable :: SigF(:)
      real(dp), allocatable :: NuSigF(:)  ! Fission Production [cm^-1]
      real(dp), allocatable :: Chi(:)     ! Fission Spectrum [-]
      real(dp), allocatable :: Src(:)     ! Fixed Source [cm^-2]
      real(dp), allocatable :: SigmaS(:,:)! Scattering Matrix [from, to]
  end type t_material

contains

  subroutine InitialiseMaterials(mats, n_groups, elem_mat_map, filename, printout)
      type(t_material), allocatable, intent(out) :: mats(:)
      integer, intent(in) :: n_groups
      integer, intent(in), allocatable :: elem_mat_map(:)
      character(len=*), intent(in) :: filename
      logical, intent(in) :: printout
      integer :: i, id, max_id, u, ios
      integer, allocatable :: all_ids(:)
      character(len=1024) :: line

      open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
          print *, 'ERROR: Could not open material file: ', trim(filename)
          stop
      end if

      max_id = 0
      do
          read(u, '(A)', iostat=ios) line
          if (ios /= 0) exit
          line = adjustl(trim(line))

          if (line == '' .or. line(1:1) == '-' .or. index(line, 'matID') > 0) cycle

          read(line, *, iostat=ios) id
          if (ios == 0) max_id = max(max_id, id)
      end do
      close(u)

      if (max_id == 0) then
          print *, 'ERROR: No materials found in file: ', trim(filename)
          stop
      end if

      allocate(mats(max_id))
      allocate(all_ids(max_id))
      all_ids = [(i, i=1, max_id)]

      call ParseMaterialDeck(mats, all_ids, n_groups, filename)

      do i = 1, max_id
          if (allocated(mats(i)%SigA)) call UpdateComputables(mats(i), n_groups)
      end do

      if (printout) call PrintMaterialSummary(mats, elem_mat_map, n_groups)
      deallocate(all_ids)
  end subroutine InitialiseMaterials

  subroutine ParseMaterialDeck(mats, requested_ids, n_groups, filename)
      type(t_material), intent(inout) :: mats(:)
      integer, intent(in)             :: requested_ids(:), n_groups
      character(len=*), intent(in)    :: filename
      integer :: u, id, g, ios, i, current_id
      character(len=1024) :: line
      
      open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      
      do
          read(u, '(A)', iostat=ios) line
          if (ios /= 0) exit
          line = adjustl(line)
          
          if (line == '' .or. line(1:1) == '-' .or. index(line, 'matID') > 0) cycle

          read(line, *, iostat=ios) current_id
          if (ios /= 0) cycle
          
          if (current_id <= 0 .or. current_id > size(mats)) cycle
          if (.not. any(requested_ids == current_id)) cycle

          if (.not. allocated(mats(current_id)%SigA)) then
              allocate(mats(current_id)%SigA(n_groups), mats(current_id)%Nu(n_groups),    &
                       mats(current_id)%SigF(n_groups), mats(current_id)%NuSigF(n_groups),&
                       mats(current_id)%Chi(n_groups),  mats(current_id)%Src(n_groups),   &
                       mats(current_id)%SigmaS(n_groups, n_groups),               &
                       mats(current_id)%D(n_groups), mats(current_id)%SigmaT(n_groups),   &
                       mats(current_id)%SigmaR(n_groups))
              mats(current_id)%SigmaS = 0.0_dp
              mats(current_id)%name   = 'Unknown'
          end if

          do g = 1, n_groups
              if (g > 1) then
                  read(u, '(A)', iostat=ios) line
                  if (ios /= 0) exit
              end if
              
              if (g == 1) then
                  read(line, *, iostat=ios) id, mats(current_id)%SigA(g), mats(current_id)%Nu(g), &
                                mats(current_id)%SigF(g), mats(current_id)%Chi(g),    &
                                mats(current_id)%Src(g), (mats(current_id)%SigmaS(g, i), i=1, n_groups), &
                                mats(current_id)%name
              else
                  read(line, *, iostat=ios) id, mats(current_id)%SigA(g), mats(current_id)%Nu(g), &
                                mats(current_id)%SigF(g), mats(current_id)%Chi(g),    &
                                mats(current_id)%Src(g), (mats(current_id)%SigmaS(g, i), i=1, n_groups)
              end if
          end do
      end do
      close(u)
  end subroutine ParseMaterialDeck

  subroutine UpdateComputables(mat, n_groups)
      type(t_material), intent(inout) :: mat
      integer, intent(in) :: n_groups
      integer :: g

      if (sum(mat%Chi) > 1e-12_dp) mat%Chi = mat%Chi / sum(mat%Chi)

      do g = 1, n_groups
          mat%NuSigF(g) = mat%Nu(g) * mat%SigF(g)
          mat%SigmaT(g) = mat%SigA(g) + sum(mat%SigmaS(g, :))
          mat%SigmaR(g) = mat%SigmaT(g) - mat%SigmaS(g, g)
          
          if (mat%SigmaT(g) > 1e-10_dp) then
              mat%D(g) = 1.0_dp / (3.0_dp * mat%SigmaT(g))
          else
              mat%D(g) = 1.0_dp / (3.0_dp * 1e-10_dp) ! Avoid division by zero
          end if
      end do
  end subroutine UpdateComputables

subroutine PrintMaterialSummary(mats, elem_mat_map, n_groups)
      type(t_material), intent(in) :: mats(:) 
      integer, intent(in), optional :: elem_mat_map(:) ! Made optional for flexibility
      integer, intent(in)          :: n_groups
      integer                      :: i, g
      character(len=64)            :: fmt_line, fmt_scat
      logical                      :: is_used

      write(fmt_scat, '(A, I2, A)') '(', n_groups, 'ES11.2)'

      print *, ""
      print '(A)', " [ MATS   ] :: ACTIVE MATERIAL DATABASE"
      print '(A)', " ========================================================================================================================================================================="
      print '(A4,2X,A15,3X,A3,4X,A10,2X,A10,2X,A10,2X,A10,2X,A10,3X,A,A)', &
            "ID", "Name", "Grp", "  SigA  ", "   Nu   ", "  SigF  ", " NuSigF ", " SigmaR ", "|", " Scattering (to g=1..N)"
      print '(A)', " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

      do i = 1, size(mats)
          if (.not. allocated(mats(i)%SigA)) cycle
          
          is_used = .true.
          if (present(elem_mat_map)) then
              if (.not. any(elem_mat_map == i)) is_used = .false.
          end if

          if (is_used) then
              do g = 1, n_groups
                  if (g == 1) then
                      write(fmt_line, '(A, A, A)') '(I4, 2X, A15, 4X, I3, 1X, 5ES12.4, 3X, A, ', trim(fmt_scat), ')'
                      print fmt_line, &
                            i, mats(i)%name, g, &
                            mats(i)%SigA(g), mats(i)%Nu(g), mats(i)%SigF(g), &
                            mats(i)%NuSigF(g), mats(i)%SigmaR(g), "|", mats(i)%SigmaS(g, :)
                  else
                      write(fmt_line, '(A, A, A)') '(25X, I3, 1X, 5ES12.4, 3X, A, ', trim(fmt_scat), ')'
                      print fmt_line, &
                            g, &
                            mats(i)%SigA(g), mats(i)%Nu(g), mats(i)%SigF(g), &
                            mats(i)%NuSigF(g), mats(i)%SigmaR(g), "|", mats(i)%SigmaS(g, :)
                  end if
              end do
              print '(A)', " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
          end if
      end do
      print *, ""
      
  end subroutine PrintMaterialSummary
end module m_materials
