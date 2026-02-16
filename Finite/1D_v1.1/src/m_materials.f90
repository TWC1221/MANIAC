module m_materials
  use m_constants
  implicit none
  private
  public :: t_material
  public :: init_material_const
  public :: SigmaA_of, SigmaS_of, SigmanuF_of, chi_of, SigmaT_of, SigmaR_of, D_of, geo_weight

  type :: t_material
    integer :: G         = 0
    integer :: geom_m    = 0      ! 0 slab, 1 cyl, 2 sph
    real(dp) :: tiny_r   = 1.0e-14_dp
    real(dp), allocatable :: SigmaA(:)      ! [G]
    real(dp), allocatable :: SigmanuF(:)    ! [G]
    real(dp), allocatable :: chi(:)         ! [G]
    real(dp), allocatable :: SigmaS(:,:)    ! [G,G]
  contains
    procedure :: SigmaA_of
    procedure :: SigmaS_of
    procedure :: SigmanuF_of
    procedure :: chi_of
    procedure :: SigmaT_of
    procedure :: SigmaR_of
    procedure :: D_of
    procedure :: geo_weight
  end type t_material

contains

  subroutine init_material_const(mats, G, geom_m)
    type(t_material), intent(out) :: mats
    integer,          intent(in)  :: G, geom_m

    mats%G      = G
    mats%geom_m = geom_m
    allocate(mats%SigmaA(G), mats%SigmanuF(G), mats%chi(G), mats%SigmaS(G,G))

    mats%SigmaA   = [0.015_dp, 0.04_dp, 0.12_dp]
    mats%SigmanuF = [0.02_dp, 0.10_dp, 0.35_dp]
    mats%chi      = [1.0_dp, 0.0_dp, 0.0_dp]
    mats%SigmaS   = transpose(reshape([ &
      0.20_dp, 0.05_dp, 0.10_dp, &
      0.10_dp, 0.25_dp, 0.07_dp, &
      0.17_dp, 0.10_dp, 0.30_dp ], [G,G]))
  end subroutine init_material_const

  pure real(dp) function SigmaA_of(this, x, g) result(Sa)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    integer,           intent(in) :: g
    Sa = this%SigmaA(g)
  end function SigmaA_of

  pure real(dp) function SigmaS_of(this, x, g, gp) result(Ss)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    integer,           intent(in) :: g, gp
    Ss = this%SigmaS(g,gp)
  end function SigmaS_of

  pure real(dp) function SigmanuF_of(this, x, g) result(Sf)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    integer,           intent(in) :: g
    Sf = this%SigmanuF(g)
  end function SigmanuF_of

  pure real(dp) function chi_of(this, x, g) result(c)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    integer,           intent(in) :: g
    c = this%chi(g)
  end function chi_of

  pure real(dp) function SigmaT_of(this, x, g) result(St)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    integer,           intent(in) :: g
    integer :: gp
    St = this%SigmaA_of(x,g)
    do gp = 1, this%G
      St = St + this%SigmaS_of(x,g,gp)
    end do
  end function SigmaT_of

  pure real(dp) function SigmaR_of(this, x, g) result(Sr)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    integer,           intent(in) :: g
    Sr = this%SigmaT_of(x,g) - this%SigmaS_of(x,g,g)
  end function SigmaR_of

  pure real(dp) function D_of(this, x, g) result(D)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    integer,           intent(in) :: g
    real(dp), parameter :: TINY = 1.0e-12_dp
    real(dp) :: St
    St = this%SigmaT_of(x,g)
    D  = 1.0_dp / (3.0_dp * max(St, TINY))
  end function D_of

  pure real(dp) function geo_weight(this, x) result(w)
    class(t_material), intent(in) :: this
    real(dp),          intent(in) :: x
    w = max(x, this%tiny_r)**this%geom_m
  end function geo_weight

end module m_materials
