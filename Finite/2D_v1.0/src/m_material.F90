module m_material
  use m_constants
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
  subroutine InitialiseMaterials(mats, n_groups, printout)
      type(t_material), allocatable, intent(out) :: mats(:)
      integer, intent(in) :: n_groups
      logical, intent(in) :: printout
      integer :: g, i, k
      real(dp) :: ref_SigA(3), ref_NuSigF(3), ref_Chi(3), ref_SigS(3,3)

      allocate(mats(1))

      ref_SigA   = [0.008_dp, 0.050_dp, 0.350_dp]
      ref_NuSigF = [0.006_dp, 0.045_dp, 0.850_dp]
      ref_Chi    = [0.970_dp, 0.030_dp, 0.000_dp]
      
      ! ref_SigS(to_group, from_group)
      ref_SigS = reshape([ &
          0.1000_dp, 0.0500_dp, 0.0100_dp, & ! Col 1: From G1 (Downscatters to 2 and 3)
          0.0020_dp, 0.2000_dp, 0.0800_dp, & ! Col 2: From G2 (Mostly 2->2 and 2->3)
          0.0000_dp, 0.100_dp, 0.8000_dp  & ! Col 3: From G3 (Mostly 3->3, tiny upscatter to 2)
          ], [3, 3])

      mats(1)%name = "Fuel"
      allocate(mats(1)%D(n_groups), mats(1)%SigmaT(n_groups), &
               mats(1)%SigmaR(n_groups), mats(1)%NuSigF(n_groups), &
               mats(1)%Chi(n_groups), mats(1)%SigmaS(n_groups, n_groups))

      mats(1)%NuSigF = ref_NuSigF(1:n_groups)
      mats(1)%Chi    = ref_Chi(1:n_groups)
      
      if (sum(mats(1)%Chi) > 0.0_dp) mats(1)%Chi = mats(1)%Chi / sum(mats(1)%Chi)

      do g = 1, n_groups
          mats(1)%SigmaS(:, g) = ref_SigS(1:n_groups, g)
          mats(1)%SigmaT(g) = ref_SigA(g) + sum(ref_SigS(:, g))
          mats(1)%SigmaR(g) = mats(1)%SigmaT(g) - mats(1)%SigmaS(g, g)
          mats(1)%D(g) = 1.0_dp / (3.0_dp * mats(1)%SigmaT(g))
      end do

      ! mats(2)%name = "Moderator"
      ! allocate(mats(2)%D(n_groups), mats(2)%SigmaT(n_groups), &
      !          mats(2)%SigmaR(n_groups), mats(2)%NuSigF(n_groups), &
      !          mats(2)%Chi(n_groups), mats(2)%SigmaS(n_groups, n_groups))

      ! mats(2)%SigA   = [0.0005_dp, 0.001_dp, 0.020_dp]
      ! mats(2)%NuSigF = [0.0_dp, 0.0_dp, 0.0_dp]
      ! mats(2)%Chi    = [0.0_dp, 0.0_dp, 0.0_dp]
      
      ! mats(2)%SigmaS = reshape([ &
      !     0.080_dp, 0.000_dp, 0.000_dp, & 
      !     0.150_dp, 0.450_dp, 0.012_dp, & 
      !     0.000_dp, 0.280_dp, 2.100_dp  & 
      !     ], [n_groups, n_groups], order=[2, 1])

      ! do g = 1, n_groups
      !     mats(2)%SigmaT(g) = mats(2)%SigA(g) + sum(mats(2)%SigmaS(:, g))
      !     mats(2)%SigmaR(g) = mats(2)%SigmaT(g) - mats(2)%SigmaS(g, g)
      !     mats(2)%D(g)      = 1.0_dp / (3.0_dp * mats(2)%SigmaT(g))
      ! end do

      if (printout) then
        write(*, '(/, A)') "=========================================================================="
        write(*, '(A, I2, A)') " MATERIAL DATA SUMMARY (", n_groups, " Energy Groups)"
        write(*, '(A)') "=========================================================================="
        
        do i = 1, size(mats)
            write(*, '(/, A, A)') " Material ID: ", mats(i)%name
            write(*, '(A)') " --------------------------------------------------------------------------"
            write(*, '(A, 10I12)') " Group:          ", (g, g=1, n_groups)
            write(*, '(A, 10F12.5)') " Diff. Coeff (D):", (mats(i)%D(g), g=1, n_groups)
            write(*, '(A, 10F12.5)') " Sigma Total:    ", (mats(i)%SigmaT(g), g=1, n_groups)
            write(*, '(A, 10F12.5)') " Sigma Removal:  ", (mats(i)%SigmaR(g), g=1, n_groups)
            write(*, '(A, 10F12.5)') " Sigma Absorbtion:", (ref_SigA(g), g=1, n_groups)

            write(*, '(A, 10F12.5)') " Nu * Sigma F:   ", (mats(i)%NuSigF(g), g=1, n_groups)
            write(*, '(A, 10F12.5)') " Fission Spect.: ", (mats(i)%Chi(g), g=1, n_groups)
            
            write(*, '(A)') " Scattering Matrix (Sigma S [g -> g']):"
            do g = 1, n_groups
                write(*, '(A, I2, A, 10F12.5)') "   from G", g, " ->", (mats(i)%SigmaS(k, g), k=1, n_groups)
            end do
        end do
        write(*, '(/, A)') "=========================================================================="
      end if

  end subroutine InitialiseMaterials
end module m_material