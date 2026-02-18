module m_material
  use m_constants
  implicit none

  type t_material
      character(len=32) :: name
      real(dp), allocatable :: D(:)       ! Diffusion coefficient [cm]
      real(dp), allocatable :: SigmaT(:)  ! Total Cross Section [cm^-1]
      real(dp), allocatable :: SigmaR(:)  ! Removal Cross Section [cm^-1]
      real(dp), allocatable :: NuSigF(:)  ! Fission Production [cm^-1]
      real(dp), allocatable :: Chi(:)     ! Fission Spectrum [-]
      real(dp), allocatable :: SigmaS(:,:)! Scattering Matrix [to, from]
  end type t_material

  contains
  subroutine InitialiseMaterials(mats, n_groups)
      type(t_material), allocatable, intent(out) :: mats(:)
      integer, intent(in) :: n_groups
      real(dp) :: SigA(n_groups)
      integer :: g

      allocate(mats(2))

      mats(1)%name = "Fuel"
      allocate(mats(1)%D(n_groups), mats(1)%SigmaT(n_groups), &
               mats(1)%SigmaR(n_groups), mats(1)%NuSigF(n_groups), &
               mats(1)%Chi(n_groups), mats(1)%SigmaS(n_groups, n_groups))

      SigA           = [0.012_dp, 0.085_dp, 0.150_dp]
      mats(1)%NuSigF = [0.006_dp, 0.045_dp, 0.210_dp]
      mats(1)%Chi    = [0.970_dp, 0.030_dp, 0.000_dp]
      
      mats(1)%SigmaS = reshape([ &
          0.1150_dp, 0.0000_dp, 0.0000_dp, & ! To G1 (from 1, 2, 3)
          0.0450_dp, 0.3500_dp, 0.0005_dp, & ! To G2 (from 1, 2, 3)
          0.0000_dp, 0.0250_dp, 0.6100_dp  & ! To G3 (from 1, 2, 3)
          ], [n_groups, n_groups], order=[2, 1])

      do g = 1, n_groups
          mats(1)%SigmaT(g) = SigA(g) + sum(mats(1)%SigmaS(:, g))
          mats(1)%SigmaR(g) = mats(1)%SigmaT(g) - mats(1)%SigmaS(g, g)
          mats(1)%D(g)      = 1.0_dp / (3.0_dp * mats(1)%SigmaT(g))
      end do

      mats(2)%name = "Moderator"
      allocate(mats(2)%D(n_groups), mats(2)%SigmaT(n_groups), &
               mats(2)%SigmaR(n_groups), mats(2)%NuSigF(n_groups), &
               mats(2)%Chi(n_groups), mats(2)%SigmaS(n_groups, n_groups))

      SigA           = [0.0005_dp, 0.001_dp, 0.020_dp]
      mats(2)%NuSigF = [0.0_dp, 0.0_dp, 0.0_dp]
      mats(2)%Chi    = [0.0_dp, 0.0_dp, 0.0_dp]
      
      mats(2)%SigmaS = reshape([ &
          0.080_dp, 0.000_dp, 0.000_dp, & 
          0.150_dp, 0.450_dp, 0.012_dp, & 
          0.000_dp, 0.280_dp, 2.100_dp  & 
          ], [n_groups, n_groups], order=[2, 1])

      do g = 1, n_groups
          mats(2)%SigmaT(g) = SigA(g) + sum(mats(2)%SigmaS(:, g))
          mats(2)%SigmaR(g) = mats(2)%SigmaT(g) - mats(2)%SigmaS(g, g)
          mats(2)%D(g)      = 1.0_dp / (3.0_dp * mats(2)%SigmaT(g))
      end do

  end subroutine InitialiseMaterials
end module m_material