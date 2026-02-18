module m_material
  use m_constants
  implicit none

  type t_material
      character(len=32) :: name
      real(dp), allocatable :: D(:)       ! Diffusion coefficient [cm]
      real(dp), allocatable :: SigmaR(:)  ! Removal Cross Section [cm^-1]
      real(dp), allocatable :: NuSigF(:)  ! Fission Production [cm^-1]
  end type t_material

  contains
  subroutine InitialiseMaterials(mats, n_groups)
      type(t_material), allocatable, intent(out) :: mats(:)
      integer, intent(in) :: n_groups
      
      allocate(mats(3))

      ! Material ID 1: Moderator
      mats(1)%name = "Moderator"
      allocate(mats(1)%D(n_groups), mats(1)%SigmaR(n_groups), mats(1)%NuSigF(n_groups))
      mats(1)%D(1)      = 0.8333_dp
      mats(1)%SigmaR(1) = 0.4_dp
      mats(1)%NuSigF(1) = 0.45_dp

      ! Material ID 2: Fuel
      mats(2)%name = "Fuel"
      allocate(mats(2)%D(n_groups), mats(2)%SigmaR(n_groups), mats(2)%NuSigF(n_groups))
      mats(2)%D(1)      = 1.20_dp
      mats(2)%SigmaR(1) = 0.15_dp
      mats(2)%NuSigF(1) = 0.0_dp

      ! Material ID 3: Moderator
      ! mats(3)%name = "Moderator"
      ! allocate(mats(3)%D(n_groups), mats(3)%SigmaR(n_groups), mats(3)%NuSigF(n_groups))
      ! mats(3)%D(1)      = 1.4_dp
      ! mats(3)%SigmaR(1) = 0.02_dp
      ! mats(3)%NuSigF(1) = 0.00_dp
  end subroutine InitialiseMaterials

end module m_material
