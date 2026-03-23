module m_material
  use m_constants
  use m_types
  implicit none

  type t_material
      character(len=32) :: name
      real(dp), allocatable :: D(:)       ! Diffusion coefficient [cm]
      real(dp), allocatable :: SigmaT(:)  ! Total Cross Section [cm^-1]
      real(dp), allocatable :: SigmaR(:)  ! Removal Cross Section [cm^-1]
      real(dp), allocatable :: SigA(:)    ! Absorption Cross Section [cm^-1]
      real(dp), allocatable :: NuSigF(:)  ! Fission Production [cm^-1]
      real(dp), allocatable :: Chi(:)     ! Fission Spectrum [-]
      real(dp), allocatable :: SigmaS(:,:)! Scattering Matrix [to, from]
  end type t_material

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
      integer, intent(in)             :: requested_ids(:) ! IDs found in mesh
      integer, intent(in)             :: n_groups
      character(len=*), intent(in)    :: filename
      
      integer :: unit_in, ios, g, current_mat
      character(len=256) :: line
      character(len=32)  :: keyword
      logical :: is_requested

      open(newunit=unit_in, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) stop "ERROR: Could not open material deck file!"

      current_mat = 0
      is_requested = .false.

      do
          read(unit_in, '(A)', iostat=ios) line
          if (ios /= 0) exit 

          if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
          read(line, *, iostat=ios) keyword

          select case (trim(keyword))
          case ('MATERIAL')
              read(line, *) keyword, current_mat
              
              is_requested = any(requested_ids == current_mat)
              
              if (is_requested .and. current_mat <= size(mats)) then
                  read(line, *) keyword, current_mat, mats(current_mat)%name
                  if (.not. allocated(mats(current_mat)%SigA)) then
                      allocate(mats(current_mat)%SigA(n_groups))
                      allocate(mats(current_mat)%NuSigF(n_groups))
                      allocate(mats(current_mat)%Chi(n_groups))
                      allocate(mats(current_mat)%SigmaS(n_groups, n_groups))
                      mats(current_mat)%SigmaS = 0.0_dp
                  end if
              else
                  is_requested = .false. ! Ignore subsequent lines for this mat
              end if

          case ('SIG_A', 'NUSIGF', 'CHI', 'SCAT_G1', 'SCAT_G2', 'SCAT_G3')
              if (is_requested) then
                  select case (trim(keyword))
                      case ('SIG_A');  read(line, *) keyword, mats(current_mat)%SigA(1:n_groups)
                      case ('NUSIGF'); read(line, *) keyword, mats(current_mat)%NuSigF(1:n_groups)
                      case ('CHI');    read(line, *) keyword, mats(current_mat)%Chi(1:n_groups)
                      case default;   
                          read(keyword(7:7), *) g
                          if (g <= n_groups) read(line, *) keyword, mats(current_mat)%SigmaS(1:n_groups, g)
                  end select
              end if
          end select
      end do
      close(unit_in)
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
          mat%SigmaT(g) = mat%SigA(g) + sum(mat%SigmaS(:, g))
          mat%SigmaR(g) = mat%SigmaT(g) - mat%SigmaS(g, g)
          
          if (mat%SigmaT(g) > 0.0_dp) then
              mat%D(g) = 1.0_dp / (3.0_dp * mat%SigmaT(g))
          else
              mat%D(g) = 1.0_dp ! Fallback for vacuum
          end if
      end do
  end subroutine UpdateComputables

  subroutine PrintMaterialSummary(mats, n_groups)
      type(t_material), intent(in) :: mats(:)
      integer, intent(in) :: n_groups
      integer :: i, g

      print*, " "
      print*, "================ MATERIAL SUMMARY ================"
      print '(A10, A15, A10, A10, A10)', "ID", "Name", "Group", "SigA", "NuSigF"
      
      do i = 1, size(mats)
          if (allocated(mats(i)%SigA)) then
              do g = 1, n_groups
                  if (g == 1) then
                      print '(I10, A15, I10, F10.5, F10.5)', i, mats(i)%name, g, &
                             mats(i)%SigA(g), mats(i)%NuSigF(g)
                  else
                      print '(A10, A15, I10, F10.5, F10.5)', "", "", g, &
                             mats(i)%SigA(g), mats(i)%NuSigF(g)
                  end if
              end do
              print*, "--------------------------------------------------"
          end if
      end do
  end subroutine PrintMaterialSummary

end module m_material

! allocate(mats(1))

!       ref_SigA   = [0.008_dp, 0.050_dp, 0.350_dp]
!       ref_NuSigF = [0.006_dp, 0.045_dp, 0.850_dp]
!       ref_Chi    = [0.970_dp, 0.030_dp, 0.000_dp]
      
!       ! ref_SigS(to_group, from_group)
!       ref_SigS = reshape([ &
!           0.1000_dp, 0.0500_dp, 0.0100_dp, & ! Col 1: From G1 (Downscatters to 2 and 3)
!           0.0020_dp, 0.2000_dp, 0.0800_dp, & ! Col 2: From G2 (Mostly 2->2 and 2->3)
!           0.0000_dp, 0.100_dp, 0.8000_dp  & ! Col 3: From G3 (Mostly 3->3, tiny upscatter to 2)
!           ], [3, 3])

!       mats(1)%name = "Fuel"
!       allocate(mats(1)%D(n_groups), mats(1)%SigmaT(n_groups), &
!                mats(1)%SigmaR(n_groups), mats(1)%NuSigF(n_groups), &
!                mats(1)%Chi(n_groups), mats(1)%SigmaS(n_groups, n_groups))

!       mats(1)%NuSigF = ref_NuSigF(1:n_groups)
!       mats(1)%Chi    = ref_Chi(1:n_groups)
      
!       if (sum(mats(1)%Chi) > 0.0_dp) mats(1)%Chi = mats(1)%Chi / sum(mats(1)%Chi)

!       do g = 1, n_groups
!           mats(1)%SigmaS(:, g) = ref_SigS(1:n_groups, g)
!           mats(1)%SigmaT(g) = ref_SigA(g) + sum(ref_SigS(:, g))
!           mats(1)%SigmaR(g) = mats(1)%SigmaT(g) - mats(1)%SigmaS(g, g)
!           mats(1)%D(g) = 1.0_dp / (3.0_dp * mats(1)%SigmaT(g))
!       end do

!       ! mats(2)%name = "Moderator"
!       ! allocate(mats(2)%D(n_groups), mats(2)%SigmaT(n_groups), &
!       !          mats(2)%SigmaR(n_groups), mats(2)%NuSigF(n_groups), &
!       !          mats(2)%Chi(n_groups), mats(2)%SigmaS(n_groups, n_groups))

!       ! mats(2)%SigA   = [0.0005_dp, 0.001_dp, 0.020_dp]
!       ! mats(2)%NuSigF = [0.0_dp, 0.0_dp, 0.0_dp]
!       ! mats(2)%Chi    = [0.0_dp, 0.0_dp, 0.0_dp]
      
!       ! mats(2)%SigmaS = reshape([ &
!       !     0.080_dp, 0.000_dp, 0.000_dp, & 
!       !     0.150_dp, 0.450_dp, 0.012_dp, & 
!       !     0.000_dp, 0.280_dp, 2.100_dp  & 
!       !     ], [n_groups, n_groups], order=[2, 1])

!       ! do g = 1, n_groups
!       !     mats(2)%SigmaT(g) = mats(2)%SigA(g) + sum(mats(2)%SigmaS(:, g))
!       !     mats(2)%SigmaR(g) = mats(2)%SigmaT(g) - mats(2)%SigmaS(g, g)
!       !     mats(2)%D(g)      = 1.0_dp / (3.0_dp * mats(2)%SigmaT(g))
!       ! end do
