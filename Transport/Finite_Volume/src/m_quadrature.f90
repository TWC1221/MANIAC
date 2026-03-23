module m_quadrature

!-----------------------------------------------------------------------!
!! Purpose:                                                            -!
!  Contains the subroutines for Gauss quadrature including:            -!
!    - Gauss quadrature weights and nodes                              -!
!    - Line elements                                                   -!
!    - Rectuangular elemenets                                          -!
!    - Triangular elements                                             -!
!    - Tetrahedral elements                                            -!
!    - Hexahedral elements                                             -!
!    - Full rotational symmetry level-symmetric quadrature 3D          -!
!    - RZ discretization quadratures                                   -!      
!                                                                      -!
!! Record of revisions:                                                -!
!   Date       Programmer     Description of change                    -!
!   ====       ==========     =====================                    -!
! 14/03/24      C. Jones         Original code                         -!
! 27/01/26      T. Charlton      3D LvL Sym quadratures                -!
! 28/01/26      T. Charlton      RZ (pp,qq) quadratures                -!
!-----------------------------------------------------------------------!

use m_constants
implicit none

type :: t_Quadrature
    real(dp), dimension(:), allocatable :: Xi, Eta, Zeta, W
    real(dp), allocatable :: dir(:,:)
    integer :: NoPoints
end type t_Quadrature

type t_sn_quadrature
    integer :: NoAngles, P, Q
    real(dp), dimension(:), allocatable   :: mu, eta, zeta
    real(dp), dimension(:,:), allocatable :: Angles
    real(dp), dimension(:), allocatable   :: w
    real(dp), dimension(:), allocatable   :: alpha1D
    real(dp), dimension(:,:), allocatable :: alpha2D

    integer, allocatable :: pp(:), qq(:)            ! maps: angle nn -> level p, order q
    integer, allocatable :: zeta_pp(:), qq_len(:)   ! per-level metadata

    real(dp), allocatable :: mu_pq(:,:)    ! (P,Q): μ at (pp,qq)
    real(dp), allocatable :: zeta_p(:)  ! (P,Q): ζ at (pp,qq)
    real(dp), allocatable :: w_pq(:,:)     ! (P,Q): weight at (pp,qq)
    integer, allocatable :: n_pq(:,:)

    contains
        procedure :: get_angle_weights
        procedure :: get_angle_vectors
        procedure :: get_number_of_angles

end type t_sn_quadrature

contains

function remap_points(x, a, b) result(x_remap)

    real(dp), intent(in)    :: x(:)
    real(dp), intent(in)    :: a, b

    real(dp), dimension(size(x)) :: x_remap

    x_remap = 0.5_dp * (b - a) * x + 0.5_dp * (b + a)

end function remap_points

function rescale_weights(w, a, b) result(w_rescale)

    real(dp), intent(in)    :: w(:)
    real(dp), intent(in)    :: a, b

    real(dp), dimension(size(w)) :: w_rescale

    w_rescale = 0.5_dp * (b - a) * w

end function rescale_weights 

function get_angle_weights(this) result (w)

    class(t_sn_quadrature), intent(in) :: this
    real(dp), dimension(this%NoAngles)  :: w

    w = this%w

end function get_angle_weights

function get_angle_vectors(this) result (angles)

    class(t_sn_quadrature), intent(in)   :: this
    real(dp), dimension(this%NoAngles,3)  :: angles

    angles = this%Angles

end function get_angle_vectors

function get_number_of_angles(this) result (no_angles)

    class(t_sn_quadrature), intent(in) :: this
    integer                             :: no_angles

    no_angles = this%NoAngles

end function get_number_of_angles

subroutine Get2DAngleQuadrature(sn_quad, SN, flag_Adjoint)
    type(t_sn_quadrature)  :: sn_quad
    integer,intent(in) :: SN
    logical, optional, intent(in) :: flag_Adjoint
    logical :: Adjoint = .false.
    integer                 :: i

    if (present(flag_Adjoint)) Adjoint = .true.

    if (SN == 2) then
        sn_quad%NoAngles = 4
        allocate(sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1), sn_quad%w(4))
        allocate(sn_quad%Angles(4,3))

        sn_quad%mu(1) = 0.5773502691896257_dp
        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        ! Calculate all angle pairs for Angles
        ! Top right quadrant
        if (Adjoint) then
            sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(2,:) = [sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(4,:) = [-sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1)]
        else 
            sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(2,:) = [-sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(4,:) = [sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(1)]
        end if

        sn_quad%w = 1.0_dp

    else if (SN == 4) then
        sn_quad%NoAngles = 12
        allocate(sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(2), sn_quad%w(12))
        allocate(sn_quad%Angles(12,3))
        
        sn_quad%mu(1) = 0.3500212_dp
        sn_quad%mu(2) = 0.8688903_dp

        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        if (Adjoint) then
            sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [-sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(4,:) = [sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(5,:) = [sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(6,:) = [sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(7,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(8,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(9,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(10,:) = [-sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(11,:) = [-sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(12,:) = [-sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(1)]
        else
            sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(4,:) = [-sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(5,:) = [-sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(6,:) = [-sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(7,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(8,:) = [-sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(9,:) = [-sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(10,:) = [sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(11,:) = [sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(12,:) = [sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(1)]
        end if
        
        sn_quad%w = 1.0_dp / 3.0_dp


    else if (SN == 6) then
        sn_quad%NoAngles = 24
        allocate(sn_quad%mu(3), sn_quad%eta(3), sn_quad%zeta(3), sn_quad%w(24))
        allocate(sn_quad%Angles(24,3))

        sn_quad%mu(1) = 0.2666355_dp
        sn_quad%mu(2) = 0.6815076_dp
        sn_quad%mu(3) = 0.9261808_dp
        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        ! Calculate all angle pairs for Angles
        if (Adjoint) then
        ! Top right quadrant
        sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(3)]
        sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(2)]
        sn_quad%Angles(3,:) = [-sn_quad%mu(3), -sn_quad%eta(1), sn_quad%zeta(1)]
        sn_quad%Angles(4,:) = [-sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(2)]
        sn_quad%Angles(5,:) = [-sn_quad%mu(2), -sn_quad%eta(2), sn_quad%zeta(1)]
        sn_quad%Angles(6,:) = [-sn_quad%mu(1), -sn_quad%eta(3), sn_quad%zeta(1)]

        ! Top left quadrant
        do i = 7, 12
            sn_quad%Angles(i,:) = sn_quad%Angles(i-6,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
        end do

        ! Bottom left quadrant
        do i = 13, 18
            sn_quad%Angles(i,:) = sn_quad%Angles(i-12,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
        end do

        ! Bottom right quadrant
        do i = 19, 24
            sn_quad%Angles(i,:) = sn_quad%Angles(i-18,:) * [1.0_dp, -1.0_dp, 1.0_dp]
        end do

        else
        sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(3)]
        sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(2)]
        sn_quad%Angles(3,:) = [sn_quad%mu(3), sn_quad%eta(1), sn_quad%zeta(1)]
        sn_quad%Angles(4,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(2)]
        sn_quad%Angles(5,:) = [sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(1)]
        sn_quad%Angles(6,:) = [sn_quad%mu(1), sn_quad%eta(3), sn_quad%zeta(1)]

        ! Top left quadrant 
        do i = 7, 12
            sn_quad%Angles(i,:) = sn_quad%Angles(i-6,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
        end do

        ! Bottom left quadrant
        do i = 13, 18
            sn_quad%Angles(i,:) = sn_quad%Angles(i-12,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
        end do

        ! Bottom right quadrant
        do i = 19, 24
            sn_quad%Angles(i,:) = sn_quad%Angles(i-18,:) * [1.0_dp, -1.0_dp, 1.0_dp]
        end do
        end if

        sn_quad%w(1:6) = (/0.1761263_dp, 0.1572071_dp, 0.1761263_dp, 0.1572071_dp, 0.1572071_dp, 0.1761263_dp/)
        sn_quad%w(7:12) = sn_quad%w(1:6)
        sn_quad%w(13:18) = sn_quad%w(1:6)
        sn_quad%w(19:24) = sn_quad%w(1:6)

    elseif (SN == 8) then
        sn_quad%NoAngles = 40
        allocate(sn_quad%mu(4), sn_quad%eta(4), sn_quad%zeta(4), sn_quad%w(40))
        allocate(sn_quad%Angles(40,3))

        sn_quad%mu(1) = 0.2182179_dp
        sn_quad%mu(2) = 0.5773503_dp
        sn_quad%mu(3) = 0.7867958_dp
        sn_quad%mu(4) = 0.9511897_dp
        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        ! Calculate all angle pairs for Angles
        if (Adjoint) then
        ! Top right quadrant
        sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(4)]
        sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(3)]
        sn_quad%Angles(3,:) = [-sn_quad%mu(3), -sn_quad%eta(1), sn_quad%zeta(2)]
        sn_quad%Angles(4,:) = [-sn_quad%mu(4), -sn_quad%eta(1), sn_quad%zeta(1)]
        sn_quad%Angles(5,:) = [-sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(3)]
        sn_quad%Angles(6,:) = [-sn_quad%mu(2), -sn_quad%eta(2), sn_quad%zeta(2)]
        sn_quad%Angles(7,:) = [-sn_quad%mu(3), -sn_quad%eta(2), sn_quad%zeta(1)]
        sn_quad%Angles(8,:) = [-sn_quad%mu(1), -sn_quad%eta(3), sn_quad%zeta(2)]
        sn_quad%Angles(9,:) = [-sn_quad%mu(2), -sn_quad%eta(3), sn_quad%zeta(1)]
        sn_quad%Angles(10,:) = [-sn_quad%mu(1), -sn_quad%eta(4), sn_quad%zeta(1)]

        ! Top left quadrant
        do i = 11, 20
            sn_quad%Angles(i,:) = sn_quad%Angles(i-10,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
        end do

        ! Bottom left quadrant
        do i = 21, 30
            sn_quad%Angles(i,:) = sn_quad%Angles(i-20,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
        end do

        ! Bottom right quadrant
        do i = 31, 40
            sn_quad%Angles(i,:) = sn_quad%Angles(i-30,:) * [1.0_dp, -1.0_dp, 1.0_dp]
        end do
        
        else
        ! Top right quadrant
        sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(4)]
        sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(3)]
        sn_quad%Angles(3,:) = [sn_quad%mu(3), sn_quad%eta(1), sn_quad%zeta(2)]
        sn_quad%Angles(4,:) = [sn_quad%mu(4), sn_quad%eta(1), sn_quad%zeta(1)]
        sn_quad%Angles(5,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(3)]
        sn_quad%Angles(6,:) = [sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(2)]
        sn_quad%Angles(7,:) = [sn_quad%mu(3), sn_quad%eta(2), sn_quad%zeta(1)]
        sn_quad%Angles(8,:) = [sn_quad%mu(1), sn_quad%eta(3), sn_quad%zeta(2)]
        sn_quad%Angles(9,:) = [sn_quad%mu(2), sn_quad%eta(3), sn_quad%zeta(1)]
        sn_quad%Angles(10,:) = [sn_quad%mu(1), sn_quad%eta(4), sn_quad%zeta(1)]

        ! Top left quadrant
        do i = 11, 20
            sn_quad%Angles(i,:) = sn_quad%Angles(i-10,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
        end do

        ! Bottom left quadrant
        do i = 21, 30
            sn_quad%Angles(i,:) = sn_quad%Angles(i-20,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
        end do

        ! Bottom right quadrant
        do i = 31, 40
            sn_quad%Angles(i,:) = sn_quad%Angles(i-30,:) * [1.0_dp, -1.0_dp, 1.0_dp]
        end do
        end if
        
        sn_quad%w(1:4) = (/0.1209877_dp, 0.0907407_dp, 0.0907407_dp, 0.1209877_dp/)
        sn_quad%w(5:7) = (/0.0907407_dp, 0.0925926_dp, 0.0907407_dp/)
        sn_quad%w(8:9) = (/0.0907407_dp, 0.0907407_dp/)
        sn_quad%w(10) = 0.1209877_dp

        sn_quad%w(11:20) = sn_quad%w(1:10)
        sn_quad%w(21:30) = sn_quad%w(1:10)
        sn_quad%w(31:40) = sn_quad%w(1:10)

    else if (SN == 16) then
        allocate(sn_quad%mu(8), sn_quad%eta(8), sn_quad%zeta(8), sn_quad%w(144))
        allocate(sn_quad%Angles(144,3))
        sn_quad%NoAngles = 144
        sn_quad%mu(1) = 0.1389568_dp
        sn_quad%mu(2) = 0.3922893_dp
        sn_quad%mu(3) = 0.5370966_dp
        sn_quad%mu(4) = 0.6504264_dp
        sn_quad%mu(5) = 0.7467506_dp
        sn_quad%mu(6) = 0.8319966_dp
        sn_quad%mu(7) = 0.9092855_dp
        sn_quad%mu(8) = 0.9805009_dp

        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        ! Calculate all angle pairs for Angles
        if (Adjoint) then
        ! Top right quadrant
        sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(8)]
        sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(7)]
        sn_quad%Angles(3,:) = [-sn_quad%mu(3), -sn_quad%eta(1), sn_quad%zeta(6)]
        sn_quad%Angles(4,:) = [-sn_quad%mu(4), -sn_quad%eta(1), sn_quad%zeta(5)]
        sn_quad%Angles(5,:) = [-sn_quad%mu(5), -sn_quad%eta(1), sn_quad%zeta(4)]
        sn_quad%Angles(6,:) = [-sn_quad%mu(6), -sn_quad%eta(1), sn_quad%zeta(3)]
        sn_quad%Angles(7,:) = [-sn_quad%mu(7), -sn_quad%eta(1), sn_quad%zeta(2)]
        sn_quad%Angles(8,:) = [-sn_quad%mu(8), -sn_quad%eta(1), sn_quad%zeta(1)]
        sn_quad%Angles(9,:) = [-sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(7)]
        sn_quad%Angles(10,:) = [-sn_quad%mu(2), -sn_quad%eta(2), sn_quad%zeta(6)]
        sn_quad%Angles(11,:) = [-sn_quad%mu(3), -sn_quad%eta(2), sn_quad%zeta(5)]
        sn_quad%Angles(12,:) = [-sn_quad%mu(4), -sn_quad%eta(2), sn_quad%zeta(4)]
        sn_quad%Angles(13,:) = [-sn_quad%mu(5), -sn_quad%eta(2), sn_quad%zeta(3)]
        sn_quad%Angles(14,:) = [-sn_quad%mu(6), -sn_quad%eta(2), sn_quad%zeta(2)]
        sn_quad%Angles(15,:) = [-sn_quad%mu(7), -sn_quad%eta(2), sn_quad%zeta(1)]
        sn_quad%Angles(16,:) = [-sn_quad%mu(1), -sn_quad%eta(3), sn_quad%zeta(6)]
        sn_quad%Angles(17,:) = [-sn_quad%mu(2), -sn_quad%eta(3), sn_quad%zeta(5)]
        sn_quad%Angles(18,:) = [-sn_quad%mu(3), -sn_quad%eta(3), sn_quad%zeta(4)]
        sn_quad%Angles(19,:) = [-sn_quad%mu(4), -sn_quad%eta(3), sn_quad%zeta(3)]
        sn_quad%Angles(20,:) = [-sn_quad%mu(5), -sn_quad%eta(3), sn_quad%zeta(2)]
        sn_quad%Angles(21,:) = [-sn_quad%mu(6), -sn_quad%eta(3), sn_quad%zeta(1)]
        sn_quad%Angles(22,:) = [-sn_quad%mu(1), -sn_quad%eta(4), sn_quad%zeta(5)]
        sn_quad%Angles(23,:) = [-sn_quad%mu(2), -sn_quad%eta(4), sn_quad%zeta(4)]
        sn_quad%Angles(24,:) = [-sn_quad%mu(3), -sn_quad%eta(4), sn_quad%zeta(3)]
        sn_quad%Angles(25,:) = [-sn_quad%mu(4), -sn_quad%eta(4), sn_quad%zeta(2)]
        sn_quad%Angles(26,:) = [-sn_quad%mu(5), -sn_quad%eta(4), sn_quad%zeta(1)]
        sn_quad%Angles(27,:) = [-sn_quad%mu(1), -sn_quad%eta(5), sn_quad%zeta(4)]
        sn_quad%Angles(28,:) = [-sn_quad%mu(2), -sn_quad%eta(5), sn_quad%zeta(3)]
        sn_quad%Angles(29,:) = [-sn_quad%mu(3), -sn_quad%eta(5), sn_quad%zeta(2)]
        sn_quad%Angles(30,:) = [-sn_quad%mu(4), -sn_quad%eta(5), sn_quad%zeta(1)]
        sn_quad%Angles(31,:) = [-sn_quad%mu(1), -sn_quad%eta(6), sn_quad%zeta(3)]
        sn_quad%Angles(32,:) = [-sn_quad%mu(2), -sn_quad%eta(6), sn_quad%zeta(2)]
        sn_quad%Angles(33,:) = [-sn_quad%mu(3), -sn_quad%eta(6), sn_quad%zeta(1)]
        sn_quad%Angles(34,:) = [-sn_quad%mu(1), -sn_quad%eta(7), sn_quad%zeta(2)]
        sn_quad%Angles(35,:) = [-sn_quad%mu(2), -sn_quad%eta(7), sn_quad%zeta(1)]
        sn_quad%Angles(36,:) = [-sn_quad%mu(1), -sn_quad%eta(8), sn_quad%zeta(1)]

        ! Top left quadrant
        do i = 37, 72
            sn_quad%Angles(i,:) = sn_quad%Angles(i-36,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
        end do

        ! Bottom left quadrant
        do i = 73, 108
            sn_quad%Angles(i,:) = sn_quad%Angles(i-72,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
        end do

        ! Bottom right quadrant
        do i = 109, 144
            sn_quad%Angles(i,:) = sn_quad%Angles(i-108,:) * [1.0_dp, -1.0_dp, 1.0_dp]
        end do

        else
        ! Top right quadrant
        sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(8)]
        sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(7)]
        sn_quad%Angles(3,:) = [sn_quad%mu(3), sn_quad%eta(1), sn_quad%zeta(6)]
        sn_quad%Angles(4,:) = [sn_quad%mu(4), sn_quad%eta(1), sn_quad%zeta(5)]
        sn_quad%Angles(5,:) = [sn_quad%mu(5), sn_quad%eta(1), sn_quad%zeta(4)]
        sn_quad%Angles(6,:) = [sn_quad%mu(6), sn_quad%eta(1), sn_quad%zeta(3)]
        sn_quad%Angles(7,:) = [sn_quad%mu(7), sn_quad%eta(1), sn_quad%zeta(2)]
        sn_quad%Angles(8,:) = [sn_quad%mu(8), sn_quad%eta(1), sn_quad%zeta(1)]
        sn_quad%Angles(9,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(7)]
        sn_quad%Angles(10,:) = [sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(6)]
        sn_quad%Angles(11,:) = [sn_quad%mu(3), sn_quad%eta(2), sn_quad%zeta(5)]
        sn_quad%Angles(12,:) = [sn_quad%mu(4), sn_quad%eta(2), sn_quad%zeta(4)]
        sn_quad%Angles(13,:) = [sn_quad%mu(5), sn_quad%eta(2), sn_quad%zeta(3)]
        sn_quad%Angles(14,:) = [sn_quad%mu(6), sn_quad%eta(2), sn_quad%zeta(2)]
        sn_quad%Angles(15,:) = [sn_quad%mu(7), sn_quad%eta(2), sn_quad%zeta(1)]
        sn_quad%Angles(16,:) = [sn_quad%mu(1), sn_quad%eta(3), sn_quad%zeta(6)]
        sn_quad%Angles(17,:) = [sn_quad%mu(2), sn_quad%eta(3), sn_quad%zeta(5)]
        sn_quad%Angles(18,:) = [sn_quad%mu(3), sn_quad%eta(3), sn_quad%zeta(4)]
        sn_quad%Angles(19,:) = [sn_quad%mu(4), sn_quad%eta(3), sn_quad%zeta(3)]
        sn_quad%Angles(20,:) = [sn_quad%mu(5), sn_quad%eta(3), sn_quad%zeta(2)]
        sn_quad%Angles(21,:) = [sn_quad%mu(6), sn_quad%eta(3), sn_quad%zeta(1)]
        sn_quad%Angles(22,:) = [sn_quad%mu(1), sn_quad%eta(4), sn_quad%zeta(5)]
        sn_quad%Angles(23,:) = [sn_quad%mu(2), sn_quad%eta(4), sn_quad%zeta(4)]
        sn_quad%Angles(24,:) = [sn_quad%mu(3), sn_quad%eta(4), sn_quad%zeta(3)]
        sn_quad%Angles(25,:) = [sn_quad%mu(4), sn_quad%eta(4), sn_quad%zeta(2)]
        sn_quad%Angles(26,:) = [sn_quad%mu(5), sn_quad%eta(4), sn_quad%zeta(1)]
        sn_quad%Angles(27,:) = [sn_quad%mu(1), sn_quad%eta(5), sn_quad%zeta(4)]
        sn_quad%Angles(28,:) = [sn_quad%mu(2), sn_quad%eta(5), sn_quad%zeta(3)]
        sn_quad%Angles(29,:) = [sn_quad%mu(3), sn_quad%eta(5), sn_quad%zeta(2)]
        sn_quad%Angles(30,:) = [sn_quad%mu(4), sn_quad%eta(5), sn_quad%zeta(1)]
        sn_quad%Angles(31,:) = [sn_quad%mu(1), sn_quad%eta(6), sn_quad%zeta(3)]
        sn_quad%Angles(32,:) = [sn_quad%mu(2), sn_quad%eta(6), sn_quad%zeta(2)]
        sn_quad%Angles(33,:) = [sn_quad%mu(3), sn_quad%eta(6), sn_quad%zeta(1)]
        sn_quad%Angles(34,:) = [sn_quad%mu(1), sn_quad%eta(7), sn_quad%zeta(2)]
        sn_quad%Angles(35,:) = [sn_quad%mu(2), sn_quad%eta(7), sn_quad%zeta(1)]
        sn_quad%Angles(36,:) = [sn_quad%mu(1), sn_quad%eta(8), sn_quad%zeta(1)]

        ! Top left quadrant
        do i = 37, 72
            sn_quad%Angles(i,:) = sn_quad%Angles(i-36,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
        end do

        ! Bottom left quadrant
        do i = 73, 108
            sn_quad%Angles(i,:) = sn_quad%Angles(i-72,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
        end do

        ! Bottom right quadrant
        do i = 109, 144
            sn_quad%Angles(i,:) = sn_quad%Angles(i-108,:) * [1.0_dp, -1.0_dp, 1.0_dp]
        end do
        end if

        ! Weights
        sn_quad%w(1:8) = (/0.0489872_dp, 0.0413296_dp, 0.0212326_dp, 0.0256207_dp, 0.0256207_dp, 0.0212326_dp, 0.0413296_dp, 0.0489872_dp/)
        sn_quad%w(9:15) = (/0.0413296_dp, 0.0360486_dp, 0.0144589_dp, 0.0344958_dp, 0.0144589_dp, 0.0360486_dp, 0.0413296_dp/)
        sn_quad%w(16:21) = (/0.0212326_dp, 0.0144589_dp, 0.0085179_dp, 0.0085179_dp, 0.0144589_dp, 0.0212326_dp/)
        sn_quad%w(22:26) = (/0.0256207_dp, 0.0344958_dp, 0.0085179_dp, 0.0344958_dp, 0.0256207_dp/)
        sn_quad%w(27:30) = (/0.0256207_dp, 0.0144589_dp, 0.0144589_dp, 0.0256207_dp/)
        sn_quad%w(31:33) = (/0.0212326_dp, 0.0360486_dp, 0.0212326_dp/) 
        sn_quad%w(34:35) = (/0.0413296_dp, 0.0413296_dp/)
        sn_quad%w(36) = 0.0489872_dp

        sn_quad%w(37:72) = sn_quad%w(1:36)
        sn_quad%w(73:108) = sn_quad%w(1:36)
        sn_quad%w(109:144) = sn_quad%w(1:36)

    end if

    end subroutine Get2DAngleQuadrature

subroutine GetRZAngleQuadrature(sn_quad, SN, flag_Adjoint)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  type(t_sn_quadrature)  :: sn_quad
  integer, intent(in)    :: SN
  logical, optional, intent(in) :: flag_Adjoint
  logical :: Adjoint

  integer :: N0, N, i, j, k, pp, qq, P, Q, beg, fin, M
  real(dp), parameter :: tol = 1.0e-12_dp
  integer, allocatable :: keep(:), seg_beg(:), seg_fin(:)
  integer, allocatable :: ord(:), ring(:)
  real(dp) :: zcur

  Adjoint =.false.
  if (present(flag_Adjoint)) Adjoint = flag_Adjoint

  call Get3DAngleQuadrature(sn_quad, SN, Adjoint)

  N0 = sn_quad%NoAngles
  allocate(keep(N0)); k = 0
  do i = 1, N0
    if (sn_quad%Angles(i,2) >= -tol) then
      k = k + 1; keep(k) = i
    end if
  end do
  call compact_keep(sn_quad, keep, k)   ! now NoAngles = k
  deallocate(keep)
  N = sn_quad%NoAngles

  allocate(ord(N)); do i = 1, N; ord(i) = i; end do
  call sort_idx_by_key(ord, N, sn_quad%Angles(:,3))

  allocate(seg_beg(N), seg_fin(N)); P = 0; i = 1
  do while (i <= N)
    P = P + 1
    seg_beg(P) = i
    zcur = sn_quad%Angles(ord(i),3)
    do
      i = i + 1
      if (i > N) exit
      if (abs(sn_quad%Angles(ord(i),3) - zcur) > tol) exit
    end do
    seg_fin(P) = i - 1
  end do

  if (allocated(sn_quad%qq_len)) deallocate(sn_quad%qq_len)
  allocate(sn_quad%qq_len(P)); sn_quad%qq_len = 0
  Q = 0
  do pp = 1, P
    sn_quad%qq_len(pp) = seg_fin(pp) - seg_beg(pp) + 1
    if (sn_quad%qq_len(pp) > Q) Q = sn_quad%qq_len(pp)
  end do

  if (allocated(sn_quad%pp))       deallocate(sn_quad%pp, sn_quad%qq)
  if (allocated(sn_quad%zeta_pp))  deallocate(sn_quad%zeta_pp)
  if (allocated(sn_quad%mu_pq))    deallocate(sn_quad%mu_pq, sn_quad%zeta_p, sn_quad%w_pq)
  if (allocated(sn_quad%n_pq))     deallocate(sn_quad%n_pq)

  allocate(sn_quad%pp(N), sn_quad%qq(N))
  allocate(sn_quad%zeta_pp(P))
  allocate(sn_quad%mu_pq(P,Q), sn_quad%zeta_p(P), sn_quad%w_pq(P,Q))
  allocate(sn_quad%n_pq(P,Q))

  sn_quad%pp = 0; sn_quad%qq = 0
  sn_quad%mu_pq   = 0.0_dp
  sn_quad%zeta_p = 0.0_dp
  sn_quad%w_pq    = 0.0_dp
  sn_quad%n_pq    = 0

  do pp = 1, P
    beg = seg_beg(pp); fin = seg_fin(pp); M = fin - beg + 1
    allocate(ring(M))
    do j = 1, M
      ring(j) = ord(beg + j - 1)
    end do

    call sort_idx_by_key(ring, M, sn_quad%Angles(:,1))  ! by μ

    sn_quad%zeta_pp(pp) = sn_quad%Angles(ring(1),3)     ! representative ζ for level

    do qq = 1, M
      i = ring(qq)
      sn_quad%pp(i) = pp
      sn_quad%qq(i) = qq
      sn_quad%mu_pq(pp,qq)   = sn_quad%Angles(i,1)
      sn_quad%zeta_p(pp) = sn_quad%Angles(i,3)
      sn_quad%w_pq(pp,qq)    = sn_quad%w(i)
      sn_quad%n_pq(pp,qq)    = i
    end do
    deallocate(ring)
  end do

  sn_quad%P = P
  sn_quad%Q = Q

contains
  subroutine compact_keep(snq, keep, mnew)
    type(t_sn_quadrature), intent(inout) :: snq
    integer, intent(in) :: keep(:), mnew
    real(dp), allocatable :: A(:,:), W(:)
    integer :: t
    allocate(A(mnew,3), W(mnew))
    do t = 1, mnew
      A(t,:) = snq%Angles(keep(t),:)
      W(t)   = snq%w(keep(t))
    end do
    if (allocated(snq%Angles)) deallocate(snq%Angles)
    if (allocated(snq%w))      deallocate(snq%w)
    if (allocated(snq%mu))     deallocate(snq%mu, snq%eta, snq%zeta)
    allocate(snq%Angles(mnew,3), snq%w(mnew))
    snq%Angles = A; snq%w = W; snq%NoAngles = mnew
    allocate(snq%mu(mnew), snq%eta(mnew), snq%zeta(mnew))
    snq%mu   = snq%Angles(:,1)
    snq%eta  = snq%Angles(:,2)
    snq%zeta = snq%Angles(:,3)
    deallocate(A, W)
  end subroutine compact_keep

  subroutine sort_idx_by_key(ind, m, key)     ! insertion sort, ascending
    integer, intent(inout) :: ind(:)
    integer, intent(in)    :: m
    real(dp), intent(in)   :: key(:)
    integer :: a, b, t
    real(dp) :: kval
    if (m <= 1) return
    do a = 2, m
      t = ind(a); kval = key(t); b = a - 1
      do
        if (b < 1) exit
        if (key(ind(b)) <= kval) exit
        ind(b+1) = ind(b); b = b - 1
      end do
      ind(b+1) = t
    end do
  end subroutine sort_idx_by_key
end subroutine GetRZAngleQuadrature

subroutine Get3DAngleQuadrature(sn_quad, SN, flag_Adjoint)
    type(t_sn_quadrature)  :: sn_quad
    integer,intent(in) :: SN
    logical, optional, intent(in) :: flag_Adjoint
    logical :: Adjoint = .false.
    integer                 :: i

    if (present(flag_Adjoint)) Adjoint = .true.

    if (SN == 2) then
        sn_quad%NoAngles = 8
        allocate(sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1), sn_quad%w(8))
        allocate(sn_quad%Angles(8,3))

        sn_quad%mu(1) = 0.5773502691896257_dp
        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        if (Adjoint) then
            sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(2,:) = [sn_quad%mu(1), -sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [sn_quad%mu(1), sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(4,:) = [-sn_quad%mu(1), sn_quad%eta(1), -sn_quad%zeta(1)]

            sn_quad%Angles(5:8,:) = sn_quad%Angles(1:4,:)
            sn_quad%Angles(5:8,3) = -sn_quad%Angles(1:4,3) 
        else 
            sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(2,:) = [-sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(4,:) = [sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(1)]

            sn_quad%Angles(5:8,:) = sn_quad%Angles(1:4,:)
            sn_quad%Angles(5:8,3) = -sn_quad%Angles(1:4,3) 
        end if

        sn_quad%w = 1.0_dp

    else if (SN == 4) then
        sn_quad%NoAngles = 24
        allocate(sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(2), sn_quad%w(24))
        allocate(sn_quad%Angles(24,3))
        
        sn_quad%mu(1) = 0.3500212_dp
        sn_quad%mu(2) = 0.8688903_dp

        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        if (Adjoint) then
            sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), -sn_quad%zeta(2)]
            sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [-sn_quad%mu(1), -sn_quad%eta(2), -sn_quad%zeta(1)]

            sn_quad%Angles(4,:) = [sn_quad%mu(1), -sn_quad%eta(1), -sn_quad%zeta(2)]
            sn_quad%Angles(5,:) = [sn_quad%mu(2), -sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(6,:) = [sn_quad%mu(1), -sn_quad%eta(2), -sn_quad%zeta(1)]

            sn_quad%Angles(7,:) = [sn_quad%mu(1), sn_quad%eta(1), -sn_quad%zeta(2)]
            sn_quad%Angles(8,:) = [sn_quad%mu(2), sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(9,:) = [sn_quad%mu(1), sn_quad%eta(2), -sn_quad%zeta(1)]

            sn_quad%Angles(10,:) = [-sn_quad%mu(1), sn_quad%eta(1), -sn_quad%zeta(2)]
            sn_quad%Angles(11,:) = [-sn_quad%mu(2), sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(12,:) = [-sn_quad%mu(1), sn_quad%eta(2), -sn_quad%zeta(1)]

            sn_quad%Angles(13:24,:) = sn_quad%Angles(1:12,:)
            sn_quad%Angles(13:24,3) = -sn_quad%Angles(1:12,3) 
        else
            sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(3,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(4,:) = [-sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(5,:) = [-sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(6,:) = [-sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(7,:) = [-sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(8,:) = [-sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(9,:) = [-sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(1)]

            sn_quad%Angles(10,:) = [sn_quad%mu(1), -sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(11,:) = [sn_quad%mu(2), -sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(12,:) = [sn_quad%mu(1), -sn_quad%eta(2), sn_quad%zeta(1)]
    
            sn_quad%Angles(13:24,:) = sn_quad%Angles(1:12,:)
            sn_quad%Angles(13:24,3) = -sn_quad%Angles(1:12,3) 
        end if
        
        sn_quad%w = 1.0_dp / 3.0_dp


    else if (SN == 6) then
        sn_quad%NoAngles = 48
        allocate(sn_quad%mu(3), sn_quad%eta(3), sn_quad%zeta(3), sn_quad%w(48))
        allocate(sn_quad%Angles(48,3))

        sn_quad%mu(1) = 0.2666355_dp
        sn_quad%mu(2) = 0.6815076_dp
        sn_quad%mu(3) = 0.9261808_dp
        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        ! Calculate all angle pairs for Angles
        if (Adjoint) then
            ! Top right quadrant
            sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), -sn_quad%zeta(3)]
            sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), -sn_quad%zeta(2)]
            sn_quad%Angles(3,:) = [-sn_quad%mu(3), -sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(4,:) = [-sn_quad%mu(1), -sn_quad%eta(2), -sn_quad%zeta(2)]
            sn_quad%Angles(5,:) = [-sn_quad%mu(2), -sn_quad%eta(2), -sn_quad%zeta(1)]
            sn_quad%Angles(6,:) = [-sn_quad%mu(1), -sn_quad%eta(3), -sn_quad%zeta(1)]

            ! Top left quadrant
            do i = 7, 12
                sn_quad%Angles(i,:) = sn_quad%Angles(i-6,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
            end do

            ! Bottom left quadrant
            do i = 13, 18
                sn_quad%Angles(i,:) = sn_quad%Angles(i-12,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
            end do

            ! Bottom right quadrant
            do i = 19, 24
                sn_quad%Angles(i,:) = sn_quad%Angles(i-18,:) * [1.0_dp, -1.0_dp, 1.0_dp]
            end do

            sn_quad%Angles(25:48,:) = sn_quad%Angles(1:24,:)
            sn_quad%Angles(25:48,3) = -sn_quad%Angles(1:24,3) 

        else
            sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(3)]
            sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(3,:) = [sn_quad%mu(3), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(4,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(2)]
            sn_quad%Angles(5,:) = [sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(1)]
            sn_quad%Angles(6,:) = [sn_quad%mu(1), sn_quad%eta(3), sn_quad%zeta(1)]

            ! Top left quadrant 
            do i = 7, 12
                sn_quad%Angles(i,:) = sn_quad%Angles(i-6,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
            end do

            ! Bottom left quadrant
            do i = 13, 18
                sn_quad%Angles(i,:) = sn_quad%Angles(i-12,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
            end do

            ! Bottom right quadrant
            do i = 19, 24
                sn_quad%Angles(i,:) = sn_quad%Angles(i-18,:) * [1.0_dp, -1.0_dp, 1.0_dp]
            end do

            sn_quad%Angles(25:48,:) = sn_quad%Angles(1:24,:)
            sn_quad%Angles(25:48,3) = -sn_quad%Angles(1:24,3) 
        end if

        sn_quad%w(1:6) = (/0.1761263_dp, 0.1572071_dp, 0.1761263_dp, 0.1572071_dp, 0.1572071_dp, 0.1761263_dp/)
        sn_quad%w(7:12) = sn_quad%w(1:6)
        sn_quad%w(13:18) = sn_quad%w(1:6)
        sn_quad%w(19:24) = sn_quad%w(1:6)
        sn_quad%w(25:48) = sn_quad%w(1:24)

    elseif (SN == 8) then
        sn_quad%NoAngles = 80
        allocate(sn_quad%mu(4), sn_quad%eta(4), sn_quad%zeta(4), sn_quad%w(80))
        allocate(sn_quad%Angles(80,3))

        sn_quad%mu(1) = 0.2182179_dp
        sn_quad%mu(2) = 0.5773503_dp
        sn_quad%mu(3) = 0.7867958_dp
        sn_quad%mu(4) = 0.9511897_dp
        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        ! Calculate all angle pairs for Angles
        if (Adjoint) then
            ! Top right quadrant
            sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), -sn_quad%zeta(4)]
            sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), -sn_quad%zeta(3)]
            sn_quad%Angles(3,:) = [-sn_quad%mu(3), -sn_quad%eta(1), -sn_quad%zeta(2)]
            sn_quad%Angles(4,:) = [-sn_quad%mu(4), -sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(5,:) = [-sn_quad%mu(1), -sn_quad%eta(2), -sn_quad%zeta(3)]
            sn_quad%Angles(6,:) = [-sn_quad%mu(2), -sn_quad%eta(2), -sn_quad%zeta(2)]
            sn_quad%Angles(7,:) = [-sn_quad%mu(3), -sn_quad%eta(2), -sn_quad%zeta(1)]
            sn_quad%Angles(8,:) = [-sn_quad%mu(1), -sn_quad%eta(3), -sn_quad%zeta(2)]
            sn_quad%Angles(9,:) = [-sn_quad%mu(2), -sn_quad%eta(3), -sn_quad%zeta(1)]
            sn_quad%Angles(10,:) = [-sn_quad%mu(1), -sn_quad%eta(4), -sn_quad%zeta(1)]

            ! Top left quadrant
            do i = 11, 20
                sn_quad%Angles(i,:) = sn_quad%Angles(i-10,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
            end do

            ! Bottom left quadrant
            do i = 21, 30
                sn_quad%Angles(i,:) = sn_quad%Angles(i-20,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
            end do

            ! Bottom right quadrant
            do i = 31, 40
                sn_quad%Angles(i,:) = sn_quad%Angles(i-30,:) * [1.0_dp, -1.0_dp, 1.0_dp]
            end do

            sn_quad%Angles(41:80,:) = sn_quad%Angles(1:40,:)
            sn_quad%Angles(41:80,3) = -sn_quad%Angles(1:40,3) 
        
        else
            ! Top right quadrant
            sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(4)]
            sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(3)]
            sn_quad%Angles(3,:) = [sn_quad%mu(3), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(4,:) = [sn_quad%mu(4), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(5,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(3)]
            sn_quad%Angles(6,:) = [sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(2)]
            sn_quad%Angles(7,:) = [sn_quad%mu(3), sn_quad%eta(2), sn_quad%zeta(1)]
            sn_quad%Angles(8,:) = [sn_quad%mu(1), sn_quad%eta(3), sn_quad%zeta(2)]
            sn_quad%Angles(9,:) = [sn_quad%mu(2), sn_quad%eta(3), sn_quad%zeta(1)]
            sn_quad%Angles(10,:) = [sn_quad%mu(1), sn_quad%eta(4), sn_quad%zeta(1)]

            ! Top left quadrant
            do i = 11, 20
                sn_quad%Angles(i,:) = sn_quad%Angles(i-10,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
            end do

            ! Bottom left quadrant
            do i = 21, 30
                sn_quad%Angles(i,:) = sn_quad%Angles(i-20,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
            end do

            ! Bottom right quadrant
            do i = 31, 40
                sn_quad%Angles(i,:) = sn_quad%Angles(i-30,:) * [1.0_dp, -1.0_dp, 1.0_dp]
            end do

            sn_quad%Angles(41:80,:) = sn_quad%Angles(1:40,:)
            sn_quad%Angles(41:80,3) = -sn_quad%Angles(1:40,3) 
        end if
        
        sn_quad%w(1:4) = (/0.1209877_dp, 0.0907407_dp, 0.0907407_dp, 0.1209877_dp/)
        sn_quad%w(5:7) = (/0.0907407_dp, 0.0925926_dp, 0.0907407_dp/)
        sn_quad%w(8:9) = (/0.0907407_dp, 0.0907407_dp/)
        sn_quad%w(10) = 0.1209877_dp

        sn_quad%w(11:20) = sn_quad%w(1:10)
        sn_quad%w(21:30) = sn_quad%w(1:10)
        sn_quad%w(31:40) = sn_quad%w(1:10)

        sn_quad%w(41:80) = sn_quad%w(1:40)

    else if (SN == 16) then
        allocate(sn_quad%mu(8), sn_quad%eta(8), sn_quad%zeta(8), sn_quad%w(288))
        allocate(sn_quad%Angles(288,3))
        sn_quad%NoAngles = 288
        sn_quad%mu(1) = 0.1389568_dp
        sn_quad%mu(2) = 0.3922893_dp
        sn_quad%mu(3) = 0.5370966_dp
        sn_quad%mu(4) = 0.6504264_dp
        sn_quad%mu(5) = 0.7467506_dp
        sn_quad%mu(6) = 0.8319966_dp
        sn_quad%mu(7) = 0.9092855_dp
        sn_quad%mu(8) = 0.9805009_dp

        sn_quad%eta = sn_quad%mu
        sn_quad%zeta = sn_quad%mu

        ! Calculate all angle pairs for Angles
        if (Adjoint) then
            ! Top right quadrant
            sn_quad%Angles(1,:) = [-sn_quad%mu(1), -sn_quad%eta(1), -sn_quad%zeta(8)]
            sn_quad%Angles(2,:) = [-sn_quad%mu(2), -sn_quad%eta(1), -sn_quad%zeta(7)]
            sn_quad%Angles(3,:) = [-sn_quad%mu(3), -sn_quad%eta(1), -sn_quad%zeta(6)]
            sn_quad%Angles(4,:) = [-sn_quad%mu(4), -sn_quad%eta(1), -sn_quad%zeta(5)]
            sn_quad%Angles(5,:) = [-sn_quad%mu(5), -sn_quad%eta(1), -sn_quad%zeta(4)]
            sn_quad%Angles(6,:) = [-sn_quad%mu(6), -sn_quad%eta(1), -sn_quad%zeta(3)]
            sn_quad%Angles(7,:) = [-sn_quad%mu(7), -sn_quad%eta(1), -sn_quad%zeta(2)]
            sn_quad%Angles(8,:) = [-sn_quad%mu(8), -sn_quad%eta(1), -sn_quad%zeta(1)]
            sn_quad%Angles(9,:) = [-sn_quad%mu(1), -sn_quad%eta(2), -sn_quad%zeta(7)]
            sn_quad%Angles(10,:) = [-sn_quad%mu(2), -sn_quad%eta(2), -sn_quad%zeta(6)]
            sn_quad%Angles(11,:) = [-sn_quad%mu(3), -sn_quad%eta(2), -sn_quad%zeta(5)]
            sn_quad%Angles(12,:) = [-sn_quad%mu(4), -sn_quad%eta(2), -sn_quad%zeta(4)]
            sn_quad%Angles(13,:) = [-sn_quad%mu(5), -sn_quad%eta(2), -sn_quad%zeta(3)]
            sn_quad%Angles(14,:) = [-sn_quad%mu(6), -sn_quad%eta(2), -sn_quad%zeta(2)]
            sn_quad%Angles(15,:) = [-sn_quad%mu(7), -sn_quad%eta(2), -sn_quad%zeta(1)]
            sn_quad%Angles(16,:) = [-sn_quad%mu(1), -sn_quad%eta(3), -sn_quad%zeta(6)]
            sn_quad%Angles(17,:) = [-sn_quad%mu(2), -sn_quad%eta(3), -sn_quad%zeta(5)]
            sn_quad%Angles(18,:) = [-sn_quad%mu(3), -sn_quad%eta(3), -sn_quad%zeta(4)]
            sn_quad%Angles(19,:) = [-sn_quad%mu(4), -sn_quad%eta(3), -sn_quad%zeta(3)]
            sn_quad%Angles(20,:) = [-sn_quad%mu(5), -sn_quad%eta(3), -sn_quad%zeta(2)]
            sn_quad%Angles(21,:) = [-sn_quad%mu(6), -sn_quad%eta(3), -sn_quad%zeta(1)]
            sn_quad%Angles(22,:) = [-sn_quad%mu(1), -sn_quad%eta(4), -sn_quad%zeta(5)]
            sn_quad%Angles(23,:) = [-sn_quad%mu(2), -sn_quad%eta(4), -sn_quad%zeta(4)]
            sn_quad%Angles(24,:) = [-sn_quad%mu(3), -sn_quad%eta(4), -sn_quad%zeta(3)]
            sn_quad%Angles(25,:) = [-sn_quad%mu(4), -sn_quad%eta(4), -sn_quad%zeta(2)]
            sn_quad%Angles(26,:) = [-sn_quad%mu(5), -sn_quad%eta(4), -sn_quad%zeta(1)]
            sn_quad%Angles(27,:) = [-sn_quad%mu(1), -sn_quad%eta(5), -sn_quad%zeta(4)]
            sn_quad%Angles(28,:) = [-sn_quad%mu(2), -sn_quad%eta(5), -sn_quad%zeta(3)]
            sn_quad%Angles(29,:) = [-sn_quad%mu(3), -sn_quad%eta(5), -sn_quad%zeta(2)]
            sn_quad%Angles(30,:) = [-sn_quad%mu(4), -sn_quad%eta(5), -sn_quad%zeta(1)]
            sn_quad%Angles(31,:) = [-sn_quad%mu(1), -sn_quad%eta(6), -sn_quad%zeta(3)]
            sn_quad%Angles(32,:) = [-sn_quad%mu(2), -sn_quad%eta(6), -sn_quad%zeta(2)]
            sn_quad%Angles(33,:) = [-sn_quad%mu(3), -sn_quad%eta(6), -sn_quad%zeta(1)]
            sn_quad%Angles(34,:) = [-sn_quad%mu(1), -sn_quad%eta(7), -sn_quad%zeta(2)]
            sn_quad%Angles(35,:) = [-sn_quad%mu(2), -sn_quad%eta(7), -sn_quad%zeta(1)]
            sn_quad%Angles(36,:) = [-sn_quad%mu(1), -sn_quad%eta(8), -sn_quad%zeta(1)]

            ! Top left quadrant
            do i = 37, 72
                sn_quad%Angles(i,:) = sn_quad%Angles(i-36,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
            end do

            ! Bottom left quadrant
            do i = 73, 108
                sn_quad%Angles(i,:) = sn_quad%Angles(i-72,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
            end do

            ! Bottom right quadrant
            do i = 109, 144
                sn_quad%Angles(i,:) = sn_quad%Angles(i-108,:) * [1.0_dp, -1.0_dp, 1.0_dp]
            end do

            sn_quad%Angles(145:288,:) = sn_quad%Angles(1:144, :)
            sn_quad%Angles(145:288,3) = -sn_quad%Angles(1:144,3)

        else
            ! Top right quadrant
            sn_quad%Angles(1,:) = [sn_quad%mu(1), sn_quad%eta(1), sn_quad%zeta(8)]
            sn_quad%Angles(2,:) = [sn_quad%mu(2), sn_quad%eta(1), sn_quad%zeta(7)]
            sn_quad%Angles(3,:) = [sn_quad%mu(3), sn_quad%eta(1), sn_quad%zeta(6)]
            sn_quad%Angles(4,:) = [sn_quad%mu(4), sn_quad%eta(1), sn_quad%zeta(5)]
            sn_quad%Angles(5,:) = [sn_quad%mu(5), sn_quad%eta(1), sn_quad%zeta(4)]
            sn_quad%Angles(6,:) = [sn_quad%mu(6), sn_quad%eta(1), sn_quad%zeta(3)]
            sn_quad%Angles(7,:) = [sn_quad%mu(7), sn_quad%eta(1), sn_quad%zeta(2)]
            sn_quad%Angles(8,:) = [sn_quad%mu(8), sn_quad%eta(1), sn_quad%zeta(1)]
            sn_quad%Angles(9,:) = [sn_quad%mu(1), sn_quad%eta(2), sn_quad%zeta(7)]
            sn_quad%Angles(10,:) = [sn_quad%mu(2), sn_quad%eta(2), sn_quad%zeta(6)]
            sn_quad%Angles(11,:) = [sn_quad%mu(3), sn_quad%eta(2), sn_quad%zeta(5)]
            sn_quad%Angles(12,:) = [sn_quad%mu(4), sn_quad%eta(2), sn_quad%zeta(4)]
            sn_quad%Angles(13,:) = [sn_quad%mu(5), sn_quad%eta(2), sn_quad%zeta(3)]
            sn_quad%Angles(14,:) = [sn_quad%mu(6), sn_quad%eta(2), sn_quad%zeta(2)]
            sn_quad%Angles(15,:) = [sn_quad%mu(7), sn_quad%eta(2), sn_quad%zeta(1)]
            sn_quad%Angles(16,:) = [sn_quad%mu(1), sn_quad%eta(3), sn_quad%zeta(6)]
            sn_quad%Angles(17,:) = [sn_quad%mu(2), sn_quad%eta(3), sn_quad%zeta(5)]
            sn_quad%Angles(18,:) = [sn_quad%mu(3), sn_quad%eta(3), sn_quad%zeta(4)]
            sn_quad%Angles(19,:) = [sn_quad%mu(4), sn_quad%eta(3), sn_quad%zeta(3)]
            sn_quad%Angles(20,:) = [sn_quad%mu(5), sn_quad%eta(3), sn_quad%zeta(2)]
            sn_quad%Angles(21,:) = [sn_quad%mu(6), sn_quad%eta(3), sn_quad%zeta(1)]
            sn_quad%Angles(22,:) = [sn_quad%mu(1), sn_quad%eta(4), sn_quad%zeta(5)]
            sn_quad%Angles(23,:) = [sn_quad%mu(2), sn_quad%eta(4), sn_quad%zeta(4)]
            sn_quad%Angles(24,:) = [sn_quad%mu(3), sn_quad%eta(4), sn_quad%zeta(3)]
            sn_quad%Angles(25,:) = [sn_quad%mu(4), sn_quad%eta(4), sn_quad%zeta(2)]
            sn_quad%Angles(26,:) = [sn_quad%mu(5), sn_quad%eta(4), sn_quad%zeta(1)]
            sn_quad%Angles(27,:) = [sn_quad%mu(1), sn_quad%eta(5), sn_quad%zeta(4)]
            sn_quad%Angles(28,:) = [sn_quad%mu(2), sn_quad%eta(5), sn_quad%zeta(3)]
            sn_quad%Angles(29,:) = [sn_quad%mu(3), sn_quad%eta(5), sn_quad%zeta(2)]
            sn_quad%Angles(30,:) = [sn_quad%mu(4), sn_quad%eta(5), sn_quad%zeta(1)]
            sn_quad%Angles(31,:) = [sn_quad%mu(1), sn_quad%eta(6), sn_quad%zeta(3)]
            sn_quad%Angles(32,:) = [sn_quad%mu(2), sn_quad%eta(6), sn_quad%zeta(2)]
            sn_quad%Angles(33,:) = [sn_quad%mu(3), sn_quad%eta(6), sn_quad%zeta(1)]
            sn_quad%Angles(34,:) = [sn_quad%mu(1), sn_quad%eta(7), sn_quad%zeta(2)]
            sn_quad%Angles(35,:) = [sn_quad%mu(2), sn_quad%eta(7), sn_quad%zeta(1)]
            sn_quad%Angles(36,:) = [sn_quad%mu(1), sn_quad%eta(8), sn_quad%zeta(1)]

            ! Top left quadrant
            do i = 37, 72
                sn_quad%Angles(i,:) = sn_quad%Angles(i-36,:) * [-1.0_dp, 1.0_dp, 1.0_dp]
            end do

            ! Bottom left quadrant
            do i = 73, 108
                sn_quad%Angles(i,:) = sn_quad%Angles(i-72,:) * [-1.0_dp, -1.0_dp, 1.0_dp]
            end do

            ! Bottom right quadrant
            do i = 109, 144
                sn_quad%Angles(i,:) = sn_quad%Angles(i-108,:) * [1.0_dp, -1.0_dp, 1.0_dp]
            end do

            ! ZETA Top right quadrant
            do i = 145, 180
                sn_quad%Angles(i,:) = sn_quad%Angles(i-144,:) * [1.0_dp, 1.0_dp, -1.0_dp]
            end do

            ! ZETA Top left quadrant
            do i = 181, 216
                sn_quad%Angles(i,:) = sn_quad%Angles(i-180,:) * [-1.0_dp, 1.0_dp, -1.0_dp]
            end do

            ! ZETA Bottom left quadrant
            do i = 217, 252
                sn_quad%Angles(i,:) = sn_quad%Angles(i-216,:) * [-1.0_dp, -1.0_dp, -1.0_dp]
            end do

            ! ZETA Bottom right quadrant
            do i = 253, 288
                sn_quad%Angles(i,:) = sn_quad%Angles(i-252,:) * [1.0_dp, -1.0_dp, -1.0_dp]
            end do
            end if

            ! Weights
            sn_quad%w(1:8) = (/0.0489872_dp, 0.0413296_dp, 0.0212326_dp, 0.0256207_dp, 0.0256207_dp, 0.0212326_dp, 0.0413296_dp, 0.0489872_dp/)
            sn_quad%w(9:15) = (/0.0413296_dp, 0.0360486_dp, 0.0144589_dp, 0.0344958_dp, 0.0144589_dp, 0.0360486_dp, 0.0413296_dp/)
            sn_quad%w(16:21) = (/0.0212326_dp, 0.0144589_dp, 0.0085179_dp, 0.0085179_dp, 0.0144589_dp, 0.0212326_dp/)
            sn_quad%w(22:26) = (/0.0256207_dp, 0.0344958_dp, 0.0085179_dp, 0.0344958_dp, 0.0256207_dp/)
            sn_quad%w(27:30) = (/0.0256207_dp, 0.0144589_dp, 0.0144589_dp, 0.0256207_dp/)
            sn_quad%w(31:33) = (/0.0212326_dp, 0.0360486_dp, 0.0212326_dp/) 
            sn_quad%w(34:35) = (/0.0413296_dp, 0.0413296_dp/)
            sn_quad%w(36) = 0.0489872_dp

            sn_quad%w(37:72) = sn_quad%w(1:36)
            sn_quad%w(73:108) = sn_quad%w(1:36)
            sn_quad%w(109:144) = sn_quad%w(1:36)
            sn_quad%w(145:180) = sn_quad%w(1:36)
            sn_quad%w(181:216) = sn_quad%w(1:36)
            sn_quad%w(217:252) = sn_quad%w(1:36)
            sn_quad%w(253:288) = sn_quad%w(1:36)
    end if

    end subroutine Get3DAngleQuadrature

subroutine GetLineQuad(Quad, IntegOrder, flag_adjoint)
    type(t_Quadrature)   :: Quad
    integer, intent(in)  :: IntegOrder
    logical, intent(in), optional :: flag_adjoint

    allocate(Quad%Xi(IntegOrder+1))
    allocate(Quad%W(IntegOrder+1))
    Quad%NoPoints = IntegOrder+1

    select case (IntegOrder)
    case(0)
        Quad%Xi(1) = 0.0_dp
        Quad%W(1) = 2.0_dp

    case(1)
        Quad%Xi(1) = -1.0_dp / sqrt(3.0_dp)
        Quad%Xi(2) = 1.0_dp / sqrt(3.0_dp)

        Quad%W(1) = 1.0_dp
        Quad%W(2) = 1.0_dp

    case(2)
        Quad%Xi(1) = - sqrt(3.0_dp / 5.0_dp)
        Quad%Xi(2) = 0.0_dp
        Quad%Xi(3) = sqrt(3.0_dp / 5.0_dp)

        Quad%W(1) = 5.0_dp / 9.0_dp
        Quad%W(2) = 8.0_dp / 9.0_dp
        Quad%W(3) = 5.0_dp / 9.0_dp

    case(3)
        Quad%Xi(1) = - 0.861136311594053_dp
        Quad%Xi(2) = - 0.339981043584856_dp
        Quad%Xi(3) = 0.339981043584856_dp
        Quad%Xi(4) = 0.861136311594053_dp

        Quad%W(1) = 0.347854845137454_dp
        Quad%W(2) = 0.652145154862546_dp
        Quad%W(3) = 0.652145154862546_dp
        Quad%W(4) = 0.347854845137454_dp

    case(4)
        Quad%Xi(1) = 0.9061798459386640_dp
        Quad%Xi(2) = 0.5384693101056831_dp
        Quad%Xi(3) = 0.0_dp
        Quad%Xi(4) = - 0.5384693101056831_dp
        Quad%Xi(5) = - 0.9061798459386640_dp

        Quad%W(1) = 0.236926885056189_dp
        Quad%W(2) = 0.478628670499366_dp
        Quad%W(3) = 0.568888888888889_dp
        Quad%W(4) = 0.478628670499366_dp
        Quad%W(5) = 0.236926885056189_dp

    case(5)
        Quad%Xi(1) = 0.932469514203152_dp
        Quad%Xi(2) = 0.661209386466265_dp
        Quad%Xi(3) = 0.238619186083197_dp
        Quad%Xi(4) = - 0.238619186083197_dp
        Quad%Xi(5) = - 0.661209386466265_dp
        Quad%Xi(6) = - 0.932469514203152_dp

        Quad%W(1) = 0.171324492379170_dp
        Quad%W(2) = 0.360761573048139_dp
        Quad%W(3) = 0.467913934572691_dp
        Quad%W(4) = 0.467913934572691_dp
        Quad%W(5) = 0.360761573048139_dp
        Quad%W(6) = 0.171324492379170_dp

    case(6)
        Quad%Xi(1) = - 0.9491079123427585_dp
        Quad%Xi(2) = - 0.7415311855993945_dp
        Quad%Xi(3) = - 0.4058451513773972_dp
        Quad%Xi(4) = 0.0_dp
        Quad%Xi(5) = 0.4058451513773972_dp
        Quad%Xi(6) = 0.7415311855993945_dp
        Quad%Xi(7) = 0.9491079123427585_dp

        Quad%W(1) = 0.1294849661688697_dp
        Quad%W(2) = 0.2797053914892766_dp
        Quad%W(3) = 0.3818300505051189_dp
        Quad%W(4) = 0.4179591836734694_dp
        Quad%W(5) = 0.3818300505051189_dp
        Quad%W(6) = 0.2797053914892766_dp
        Quad%W(7) = 0.1294849661688697_dp

    case(7)
        Quad%Xi(1) = - 0.9602898564975363_dp
        Quad%Xi(2) = - 0.7966664774136267_dp
        Quad%Xi(3) = - 0.5255324099163290_dp
        Quad%Xi(4) = - 0.1834346424956498_dp
        Quad%Xi(5) = 0.1834346424956498_dp
        Quad%Xi(6) = 0.5255324099163290_dp
        Quad%Xi(7) = 0.7966664774136267_dp
        Quad%Xi(8) = 0.9602898564975363_dp

        Quad%W(1) = 0.1012285362903763_dp
        Quad%W(2) = 0.2223810344533745_dp
        Quad%W(3) = 0.3137066458778873_dp
        Quad%W(4) = 0.3626837833783620_dp
        Quad%W(5) = 0.3626837833783620_dp
        Quad%W(6) = 0.3137066458778873_dp
        Quad%W(7) = 0.2223810344533745_dp
        Quad%W(8) = 0.1012285362903763_dp

    case(8)
        Quad%Xi(1) = - 0.9681602395076261_dp
        Quad%Xi(2) = - 0.8360311073266358_dp
        Quad%Xi(3) = - 0.6133714327005904_dp
        Quad%Xi(4) = - 0.3242534234038089_dp
        Quad%Xi(5) = 0.0_dp
        Quad%Xi(6) = 0.3242534234038089_dp
        Quad%Xi(7) = 0.6133714327005904_dp
        Quad%Xi(8) = 0.8360311073266358_dp
        Quad%Xi(9) = 0.9681602395076261_dp

        Quad%W(1) = 0.0812743883615744_dp
        Quad%W(2) = 0.1806481606948574_dp
        Quad%W(3) = 0.2606106964029354_dp
        Quad%W(4) = 0.3123470770400029_dp
        Quad%W(5) = 0.3302393550012598_dp
        Quad%W(6) = 0.3123470770400029_dp
        Quad%W(7) = 0.2606106964029354_dp
        Quad%W(8) = 0.1806481606948574_dp
        Quad%W(9) = 0.0812743883615744_dp

    case(9)
        Quad%Xi(1) = - 0.9739065285171717_dp
        Quad%Xi(2) = - 0.8650633666889845_dp
        Quad%Xi(3) = - 0.6794095682990244_dp
        Quad%Xi(4) = - 0.4333953941292472_dp
        Quad%Xi(5) = - 0.1488743389816312_dp
        Quad%Xi(6) = 0.1488743389816312_dp
        Quad%Xi(7) = 0.4333953941292472_dp
        Quad%Xi(8) = 0.6794095682990244_dp
        Quad%Xi(9) = 0.8650633666889845_dp
        Quad%Xi(10) = 0.9739065285171717_dp

        Quad%W(1) = 0.0666713443086881_dp
        Quad%W(2) = 0.1494513491505806_dp
        Quad%W(3) = 0.2190863625159820_dp
        Quad%W(4) = 0.2692667193099963_dp
        Quad%W(5) = 0.2955242247147529_dp
        Quad%W(6) = 0.2955242247147529_dp
        Quad%W(7) = 0.2692667193099963_dp
        Quad%W(8) = 0.2190863625159820_dp
        Quad%W(9) = 0.1494513491505806_dp
        Quad%W(10) = 0.0666713443086881_dp

    case(10)
        Quad%Xi(1) = - 0.9782286581460570_dp
        Quad%Xi(2) = - 0.8870625997680953_dp
        Quad%Xi(3) = - 0.7301520055740494_dp
        Quad%Xi(4) = - 0.5190961292068118_dp  
        Quad%Xi(5) = - 0.2695431559523450_dp
        Quad%Xi(6) = 0.0_dp
        Quad%Xi(7) = 0.2695431559523450_dp
        Quad%Xi(8) = 0.5190961292068118_dp
        Quad%Xi(9) = 0.7301520055740494_dp
        Quad%Xi(10) = 0.8870625997680953_dp
        Quad%Xi(11) = 0.9782286581460570_dp

        Quad%W(1) = 0.0556685671161737_dp
        Quad%W(2) = 0.1255803694649046_dp
        Quad%W(3) = 0.1862902109277343_dp
        Quad%W(4) = 0.2331937645919905_dp
        Quad%W(5) = 0.2628045445102467_dp
        Quad%W(6) = 0.2729250867779006_dp
        Quad%W(7) = 0.2628045445102467_dp
        Quad%W(8) = 0.2331937645919905_dp
        Quad%W(9) = 0.1862902109277343_dp
        Quad%W(10) = 0.1255803694649046_dp
        Quad%W(11) = 0.0556685671161737_dp
    case (11)
        Quad%Xi(1) = -0.9815606342467192_dp
        Quad%Xi(2) = -0.9041172563704749_dp
        Quad%Xi(3) = -0.7699026741943047_dp
        Quad%Xi(4) = -0.5873179542866175_dp
        Quad%Xi(5) = -0.3678314989981802_dp
        Quad%Xi(6) = -0.1252334085114689_dp
        Quad%Xi(7) = 0.1252334085114689_dp
        Quad%Xi(8) = 0.3678314989981802_dp
        Quad%Xi(9) = 0.5873179542866175_dp
        Quad%Xi(10) = 0.7699026741943047_dp
        Quad%Xi(11) = 0.9041172563704749_dp
        Quad%Xi(12) = 0.9815606342467192_dp

        Quad%W(1) = 0.0471753363865118_dp
        Quad%W(2) = 0.1069393259953184_dp
        Quad%W(3) = 0.1600783285433462_dp
        Quad%W(4) = 0.2031674267230659_dp
        Quad%W(5) = 0.2334925365383548_dp
        Quad%W(6) = 0.2491470458134028_dp
        Quad%W(7) = 0.2491470458134028_dp
        Quad%W(8) = 0.2334925365383548_dp
        Quad%W(9) = 0.2031674267230659_dp
        Quad%W(10) = 0.1600783285433462_dp
        Quad%W(11) = 0.1069393259953184_dp
        Quad%W(12) = 0.0471753363865118_dp

    case(12)
        Quad%Xi(1) = -0.9841830547185881_dp
        Quad%Xi(2) = -0.9175983992229779_dp
        Quad%Xi(3) = -0.8015780907333099_dp
        Quad%Xi(4) = -0.6423493394403402_dp
        Quad%Xi(5) = -0.4484927510364469_dp
        Quad%Xi(6) = -0.2304583159551348_dp
        Quad%Xi(7) = 0.0_dp
        Quad%Xi(8) = 0.2304583159551348_dp
        Quad%Xi(9) = 0.4484927510364469_dp
        Quad%Xi(10) = 0.6423493394403402_dp
        Quad%Xi(11) = 0.8015780907333099_dp
        Quad%Xi(12) = 0.9175983992229779_dp
        Quad%Xi(13) = 0.9841830547185881_dp

        Quad%W(1) = 0.0404840047653159_dp
        Quad%W(2) = 0.0921214998377285_dp
        Quad%W(3) = 0.1388735102197872_dp
        Quad%W(4) = 0.1781459807619457_dp
        Quad%W(5) = 0.2078160475368885_dp
        Quad%W(6) = 0.2262831802628972_dp
        Quad%W(7) = 0.2325515532308739_dp
        Quad%W(8:13) = Quad%W(6:1:-1)

    case(13)
        Quad%Xi(1) = -0.9862838086968123_dp
        Quad%Xi(2) = -0.9284348836635735_dp
        Quad%Xi(3) = -0.8272013150697650_dp
        Quad%Xi(4) = -0.6872929048116855_dp
        Quad%Xi(5) = -0.5152486363581541_dp
        Quad%Xi(6) = -0.3191123689278897_dp
        Quad%Xi(7) = -0.1080549487073437_dp
        Quad%Xi(8:14) = - Quad%Xi(7:1:-1)

        Quad%W(1) = 0.0351194603317519_dp
        Quad%W(2) = 0.0801580871597602_dp
        Quad%W(3) = 0.1215185706879032_dp
        Quad%W(4) = 0.1572031671581935_dp
        Quad%W(5) = 0.1855383974779378_dp
        Quad%W(6) = 0.2051984637212956_dp
        Quad%W(7) = 0.2152638534631578_dp
        Quad%W(8:14) = Quad%W(7:1:-1)

    case(14)
        Quad%Xi(1) = -0.9879925180204854_dp
        Quad%Xi(2) = -0.9372733924007060_dp
        Quad%Xi(3) = -0.8482065834104272_dp
        Quad%Xi(4) = -0.7244177313601701_dp
        Quad%Xi(5) = -0.5709721726085388_dp
        Quad%Xi(6) = -0.3941513470775634_dp
        Quad%Xi(7) = -0.2011940939974345_dp
        Quad%Xi(8) = 0.0_dp
        Quad%Xi(9:15) = - Quad%Xi(7:1:-1)

        Quad%W(1) = 0.0307532419961173_dp
        Quad%W(2) = 0.0703660474881081_dp
        Quad%W(3) = 0.1071592204671719_dp
        Quad%W(4) = 0.1395706779261543_dp
        Quad%W(5) = 0.1662692058169939_dp
        Quad%W(6) = 0.1861610000155622_dp
        Quad%W(7) = 0.1984314853271116_dp
        Quad%W(8) = 0.2025782419255613_dp
        Quad%W(9:15) = Quad%W(7:1:-1)

    case(15)
        Quad%Xi(1) = -0.9894009349916499_dp
        Quad%Xi(2) = -0.9445750230732326_dp
        Quad%Xi(3) = -0.8656312023878318_dp
        Quad%Xi(4) = -0.7554044083550030_dp
        Quad%Xi(5) = -0.6178762444026438_dp
        Quad%Xi(6) = -0.4580167776572274_dp
        Quad%Xi(7) = -0.2816035507792589_dp
        Quad%Xi(8) = -0.0950125098376374_dp
        Quad%Xi(9:16) = - Quad%Xi(8:1:-1)

        Quad%W(1) = 0.0271524594117541_dp
        Quad%W(2) = 0.0622535239386479_dp
        Quad%W(3) = 0.0951585116824928_dp
        Quad%W(4) = 0.1246289712555339_dp
        Quad%W(5) = 0.1495959888165767_dp
        Quad%W(6) = 0.1691565193950025_dp
        Quad%W(7) = 0.1826034150449236_dp
        Quad%W(8) = 0.1894506104550685_dp
        Quad%W(9:16) = Quad%W(8:1:-1)

    case(16)
        Quad%Xi(1) = -0.9905754753144174_dp
        Quad%Xi(2) = -0.9506755217687678_dp
        Quad%Xi(3) = -0.8802391537269859_dp
        Quad%Xi(4) = -0.7815140038968014_dp
        Quad%Xi(5) = -0.6576711592166907_dp
        Quad%Xi(6) = -0.5126905370864769_dp
        Quad%Xi(7) = -0.3512317634538763_dp
        Quad%Xi(8) = -0.1784841814958478_dp
        Quad%Xi(9) = 0.0_dp
        Quad%Xi(10:17) = - Quad%Xi(8:1:-1)

        Quad%W(1) = 0.0241483028685479_dp
        Quad%W(2) = 0.0554595293739872_dp
        Quad%W(3) = 0.0850361483171792_dp
        Quad%W(4) = 0.1118838471934039_dp
        Quad%W(5) = 0.1351363684685255_dp
        Quad%W(6) = 0.1540457610768103_dp
        Quad%W(7) = 0.1680041021564500_dp
        Quad%W(8) = 0.1765627053669926_dp
        Quad%W(9) = 0.1794464703562065_dp
        Quad%W(10:17) = Quad%W(8:1:-1)

    case(17)
        Quad%Xi(1) = -0.9915651684209309_dp
        Quad%Xi(2) = -0.9558239495713977_dp
        Quad%Xi(3) = -0.8926024664975557_dp
        Quad%Xi(4) = -0.8037049589725231_dp
        Quad%Xi(5) = -0.6916870430603532_dp
        Quad%Xi(6) = -0.5597708310739475_dp
        Quad%Xi(7) = -0.4117511614628426_dp
        Quad%Xi(8) = -0.2518862256915055_dp
        Quad%Xi(9) = -0.0847750130417353_dp
        Quad%Xi(10:18) = - Quad%Xi(9:1:-1)

        Quad%W(1) = 0.0216160135264833_dp
        Quad%W(2) = 0.0497145488949698_dp
        Quad%W(3) = 0.0764257302548891_dp
        Quad%W(4) = 0.1009420441062872_dp
        Quad%W(5) = 0.1225552067114782_dp
        Quad%W(6) = 0.1406429146706507_dp
        Quad%W(7) = 0.1546846751262652_dp
        Quad%W(8) = 0.1642764837458327_dp
        Quad%W(9) = 0.1691423829631436_dp
        Quad%W(10:18) = Quad%W(9:1:-1)

    case(18)
        Quad%Xi(1) = -0.9924068438435844_dp
        Quad%Xi(2) = -0.9602081521348300_dp
        Quad%Xi(3) = -0.9031559036148179_dp
        Quad%Xi(4) = -0.8227146565371428_dp
        Quad%Xi(5) = -0.7209661773352294_dp
        Quad%Xi(6) = -0.6005453046616810_dp
        Quad%Xi(7) = -0.4645707413759600_dp
        Quad%Xi(8) = -0.3165640999636298_dp
        Quad%Xi(9) = -0.1603586456402254_dp
        Quad%Xi(10) = 0.0_dp
        Quad%Xi(11:19) = - Quad%Xi(9:1:-1)

        Quad%W(1) = 0.0194617882297265_dp
        Quad%W(2) = 0.0448142267656996_dp
        Quad%W(3) = 0.0690445427376412_dp
        Quad%W(4) = 0.0914900216224499_dp
        Quad%W(5) = 0.1115666455473339_dp
        Quad%W(6) = 0.1287539625393362_dp
        Quad%W(7) = 0.1426067021736066_dp
        Quad%W(8) = 0.1527660420658597_dp
        Quad%W(9) = 0.1589688433939543_dp
        Quad%W(10) = 0.1610544498487837_dp
        Quad%W(11:19) = Quad%W(9:1:-1)

    case(19)
        Quad%Xi(1) = -0.9931285991850949_dp
        Quad%Xi(2) = -0.9639719272779138_dp
        Quad%Xi(3) = -0.9122344282513259_dp
        Quad%Xi(4) = -0.8391169718222188_dp
        Quad%Xi(5) = -0.7463319064601508_dp
        Quad%Xi(6) = -0.6360536807265150_dp
        Quad%Xi(7) = -0.5108670019508271_dp
        Quad%Xi(8) = -0.3737060887154196_dp
        Quad%Xi(9) = -0.2277858511416451_dp
        Quad%Xi(10) = -0.0765265211334973_dp
        Quad%Xi(11:20) = - Quad%Xi(10:1:-1)

        Quad%W(1) = 0.0176140071391521_dp
        Quad%W(2) = 0.0406014298003869_dp
        Quad%W(3) = 0.0626720483341091_dp
        Quad%W(4) = 0.0832767415767048_dp
        Quad%W(5) = 0.1019301198172404_dp
        Quad%W(6) = 0.1181945319615184_dp
        Quad%W(7) = 0.1316886384491766_dp
        Quad%W(8) = 0.1420961093183820_dp
        Quad%W(9) = 0.1491729864726037_dp
        Quad%W(10) = 0.1527533871307259_dp
        Quad%W(11:20) = Quad%W(10:1:-1)

    case(20)
        Quad%Xi(1) = -0.9937521706203895_dp
        Quad%Xi(2) = -0.9672268385663063_dp
        Quad%Xi(3) = -0.9200993341504008_dp
        Quad%Xi(4) = -0.8533633645833173_dp
        Quad%Xi(5) = -0.7684399634756779_dp
        Quad%Xi(6) = -0.6671388041974123_dp
        Quad%Xi(7) = -0.5516188358872198_dp
        Quad%Xi(8) = -0.4243421202074388_dp
        Quad%Xi(9) = -0.2880213168024011_dp
        Quad%Xi(10) = -0.1455618541608951_dp
        Quad%Xi(11) = 0.0_dp
        Quad%Xi(12:21) = - Quad%Xi(10:1:-1)

        Quad%W(1) = 0.0160172282577743_dp
        Quad%W(2) = 0.0369537897708525_dp
        Quad%W(3) = 0.0571344254268572_dp
        Quad%W(4) = 0.0761001136283793_dp
        Quad%W(5) = 0.0934444234560339_dp
        Quad%W(6) = 0.1087972991671484_dp
        Quad%W(7) = 0.1218314160537285_dp
        Quad%W(8) = 0.1322689386333374_dp
        Quad%W(9) = 0.1398873947910731_dp
        Quad%W(10) = 0.1445244039899700_dp
        Quad%W(11) = 0.1460811336496904_dp
        Quad%W(12:21) = Quad%W(10:1:-1)

    end select

    if (present(flag_adjoint)) Quad%Xi = - Quad%Xi

end subroutine GetLineQuad
    
subroutine QuadrilateralQuadrature(Quad, QuadBound, IntegOrder)
    
    type(t_Quadrature)                  :: Quad, QuadBound
    integer, intent(in)                 :: IntegOrder
    integer :: ii, j, k

    allocate(Quad%Xi((IntegOrder+1)**2))
    allocate(Quad%Eta((IntegOrder+1)**2))
    allocate(Quad%W((IntegOrder+1)**2))
    Quad%NoPoints = (IntegOrder+1)**2

    k = 1
    do ii = 1, IntegOrder + 1
        do j = 1, IntegOrder + 1
            Quad%Xi(k) = QuadBound%Xi(ii)
            Quad%Eta(k) = QuadBound%Xi(j)
            Quad%W(k) = QuadBound%W(ii) * QuadBound%W(j)
            k = k + 1
        end do
    end do

end subroutine QuadrilateralQuadrature

subroutine TriangleQuadrature(Quad, IntegOrder)
    
    type(t_Quadrature)  :: Quad
    integer, intent(in) :: IntegOrder

    select case (IntegOrder)
    case(1)
        allocate(Quad%Xi(1))
        allocate(Quad%Eta(1))
        allocate(Quad%W(1))
        Quad%NoPoints = 1

        Quad%Xi(1) = 1.0_dp / 3.0_dp
        Quad%Eta(1) = 1.0_dp / 3.0_dp
        Quad%W(1) = 1.0_dp

    case(2)
        allocate(Quad%Xi(3))
        allocate(Quad%Eta(3))
        allocate(Quad%W(3))
        Quad%NoPoints = 3

        Quad%Xi(1) = 1.0_dp / 6.0_dp
        Quad%Xi(2) = 2.0_dp / 3.0_dp
        Quad%Xi(3) = 1.0_dp / 6.0_dp
    
        Quad%Eta(1) = 1.0_dp / 6.0_dp
        Quad%Eta(2) = 1.0_dp / 6.0_dp
        Quad%Eta(3) = 2.0_dp / 3.0_dp

        Quad%W(1) = 1.0_dp / 3.0_dp
        Quad%W(2) = 1.0_dp / 3.0_dp
        Quad%W(3) = 1.0_dp / 3.0_dp

    case(3)
        allocate(Quad%Xi(7))
        allocate(Quad%Eta(7))
        allocate(Quad%W(7))
        Quad%NoPoints = 7

        Quad%Xi(1) = 0.1012865073235_dp
        Quad%Xi(2) = 0.7974269853531_dp
        Quad%Xi(3) = Quad%Xi(1)
        Quad%Xi(4) = 0.4701420641051_dp
        Quad%Xi(5) = Quad%Xi(4)
        Quad%Xi(6) = 0.0597158717898_dp
        Quad%Xi(7) = 1.0_dp / 3.0_dp

        Quad%Eta(1) = Quad%Xi(1)
        Quad%Eta(2) = Quad%Xi(1)
        Quad%Eta(3) = Quad%Xi(2)
        Quad%Eta(4) = Quad%Xi(6)
        Quad%Eta(5) = Quad%Xi(4)
        Quad%Eta(6) = Quad%Xi(4)
        Quad%Eta(7) = Quad%Xi(7)

        Quad%W(1) = 0.1259391805448_dp
        Quad%W(2) = 0.1259391805448_dp
        Quad%W(3) = 0.1259391805448_dp
        Quad%W(4) = 0.1323941527885_dp
        Quad%W(5) = 0.1323941527885_dp
        Quad%W(6) = 0.1323941527885_dp
        Quad%W(7) = 0.2250000000000_dp

    case(4)
        allocate(Quad%Xi(12))
        allocate(Quad%Eta(12))
        allocate(Quad%W(12))
        Quad%NoPoints = 12

        Quad%Xi(1) = 0.063089014491502_dp
        Quad%Xi(2) = 0.873821971016996_dp
        Quad%Xi(3) = Quad%Xi(1)
        Quad%Xi(4) = 0.310352451033785_dp
        Quad%Xi(5) = 0.636502499121399_dp
        Quad%Xi(6) = 0.636502499121399_dp
        Quad%Xi(7) = Quad%Xi(4)
        Quad%Xi(8) = 0.053145049844816_dp
        Quad%Xi(9) = Quad%Xi(8)
        Quad%Xi(10) = 0.249286745170910_dp
        Quad%Xi(11) = 0.501426509658179_dp
        Quad%Xi(12) = 0.249286745170910_dp

        Quad%Eta(1) = Quad%Xi(1)
        Quad%Eta(2) = Quad%Xi(1)
        Quad%Eta(3) = Quad%Xi(2)
        Quad%Eta(4) = Quad%Xi(8)
        Quad%Eta(5) = Quad%Xi(8)
        Quad%Eta(6) = 0.310352451033785_dp
        Quad%Eta(7) = 0.636502499121399_dp
        Quad%Eta(8) = 0.636502499121399_dp
        Quad%Eta(9) = 0.310352451033785_dp
        Quad%Eta(10) = 0.249286745170910_dp
        Quad%Eta(11) = 0.249286745170910_dp
        Quad%Eta(12) = 0.501426509658179_dp

        Quad%W(1:3) = 0.050844906370207_dp
        Quad%W(4:9) = 0.082851075618374_dp
        Quad%W(10:12) = 0.116786275726379_dp

    case(5)
        allocate(Quad%Xi(16))
        allocate(Quad%Eta(16))
        allocate(Quad%W(16))
        Quad%NoPoints = 16

        Quad%Xi(1) = 1.0_dp / 3.0_dp
        Quad%Xi(2) = 0.459292588292723_dp
        Quad%Xi(3) = 0.459292588292723_dp
        Quad%Xi(4) = 0.081414823414554_dp

        Quad%Xi(5) = 0.170569307751760_dp
        Quad%Xi(6) = 0.658861384496480_dp
        Quad%Xi(7) = 0.170569307751760_dp

        Quad%Xi(8) = 0.050547228317031_dp
        Quad%Xi(9) = 0.898905543365938_dp
        Quad%Xi(10) = 0.050547228317031_dp

        Quad%Xi(11) = 0.263112829634638_dp
        Quad%Xi(12) = 0.728492392955404_dp

        Quad%Xi(13) = 0.008394777409958_dp
        Quad%Xi(14) = 0.008394777409958_dp

        Quad%Xi(15) = 0.728492392955404_dp
        Quad%Xi(16) = 0.263112829634638_dp

        Quad%Eta(1) = Quad%Xi(1)
        Quad%Eta(2) = 0.459292588292723_dp
        Quad%Eta(3) = 0.081414823414554_dp
        Quad%Eta(4) = 0.459292588292723_dp

        Quad%Eta(5) = 0.170569307751760_dp
        Quad%Eta(6) = 0.170569307751760_dp
        Quad%Eta(7) = 0.658861384496480_dp

        Quad%Eta(8) = 0.050547228317031_dp
        Quad%Eta(9) = 0.050547228317031_dp
        Quad%Eta(10) = 0.898905543365938_dp

        Quad%Eta(11) = 0.008394777409958_dp
        Quad%Eta(12) = 0.008394777409958_dp

        Quad%Eta(13) = 0.263112829634638_dp
        Quad%Eta(14) = 0.728492392955404_dp

        Quad%Eta(15) = 0.263112829634638_dp
        Quad%Eta(16) = 0.728492392955404_dp


        Quad%W(1) = 0.144315607677787_dp
        Quad%W(2:4) = 0.095091634267285_dp
        Quad%W(5:7) = 0.103217370534718_dp
        Quad%W(8:10) = 0.032458497623198_dp
        Quad%W(11:16) = 0.027230314174435_dp
    case(6)
        allocate(Quad%Xi(25))
        allocate(Quad%Eta(25))
        allocate(Quad%W(25))
        Quad%NoPoints = 25

        Quad%Xi(1) = 0.333333333333333_dp

        Quad%Xi(2) = 0.485577633383657_dp
        Quad%Xi(3) = 0.485577633383657_dp
        Quad%Xi(4) = 0.028844733232685_dp

        Quad%Xi(5) = 0.109481575485037_dp
        Quad%Xi(6) = 0.109481575485037_dp
        Quad%Xi(7) = 0.781036849029926_dp

        Quad%Xi(8) = 0.141707219414880_dp
        Quad%Xi(9) = 0.141707219414880_dp
        Quad%Xi(10) = 0.307939838764121_dp
        Quad%Xi(11) = 0.307939838764121_dp
        Quad%Xi(12) = 0.550352941820999_dp
        Quad%Xi(13) = 0.550352941820999_dp

        Quad%Xi(14) = 0.025003534762686_dp
        Quad%Xi(15) = 0.025003534762686_dp
        Quad%Xi(16) = 0.246672560639903_dp
        Quad%Xi(17) = 0.246672560639903_dp
        Quad%Xi(18) = 0.728323904597411_dp
        Quad%Xi(19) = 0.728323904597411_dp

        Quad%Xi(20) = 0.009540815400299_dp
        Quad%Xi(21) = 0.009540815400299_dp
        Quad%Xi(22) = 0.066803251012200_dp
        Quad%Xi(23) = 0.066803251012200_dp
        Quad%Xi(24) = 0.923655933587500_dp
        Quad%Xi(25) = 0.923655933587500_dp

        Quad%Eta(1) = 0.333333333333333_dp

        Quad%Eta(2) = 0.485577633383657_dp
        Quad%Eta(3) = 0.028844733232685_dp
        Quad%Eta(4) = 0.485577633383657_dp

        Quad%Eta(5) = 0.109481575485037_dp
        Quad%Eta(6) = 0.781036849029926_dp
        Quad%Eta(7) = 0.109481575485037_dp

        Quad%Eta(8) = 0.307939838764121_dp
        Quad%Eta(9) = 0.550352941820999_dp
        Quad%Eta(10) = 0.141707219414880_dp
        Quad%Eta(11) = 0.550352941820999_dp
        Quad%Eta(12) = 0.141707219414880_dp
        Quad%Eta(13) = 0.307939838764121_dp

        Quad%Eta(14) = 0.246672560639903_dp
        Quad%Eta(15) = 0.728323904597411_dp
        Quad%Eta(16) = 0.025003534762686_dp
        Quad%Eta(17) = 0.728323904597411_dp
        Quad%Eta(18) = 0.025003534762686_dp
        Quad%Eta(19) = 0.246672560639903_dp

        Quad%Eta(20) = 0.066803251012200_dp
        Quad%Eta(21) = 0.923655933587500_dp
        Quad%Eta(22) = 0.009540815400299_dp
        Quad%Eta(23) = 0.923655933587500_dp
        Quad%Eta(24) = 0.009540815400299_dp
        Quad%Eta(25) = 0.066803251012200_dp


        Quad%W(1) = 0.090817990382754_dp
        Quad%W(2:4) = 0.036725957756467_dp
        Quad%W(5:7) = 0.045321059435528_dp
        Quad%W(8:13) = 0.072757916845420_dp
        Quad%W(14:19) = 0.028327242531057_dp
        Quad%W(20:25) = 0.009421666963733_dp

    end select

    ! select case (QuadBound%NoPoints)
    !     case(1)
    !         allocate(QuadBound%Xi(1))
    !         allocate(QuadBound%Eta(1))
    !         allocate(QuadBound%W(1))

    !         QuadBound%Xi(1) = 0.5_dp
    !         QuadBound%Eta(1) = 0.5_dp
    !         QuadBound%W(1) = 1.0_dp

    !     case(2)
    !         allocate(QuadBound%Xi(2))
    !         allocate(QuadBound%W(2))

    !         QuadBound%Xi(1) = 1.0_dp / sqrt(3.0_dp)
    !         QuadBound%Xi(2) = - 1.0_dp / sqrt(3.0_dp)

    !         QuadBound%W(1) = 1.0_dp
    !         QuadBound%W(2) = 1.0_dp

    !     case(3)
    !         allocate(QuadBound%Xi(3))
    !         allocate(QuadBound%W(3))


    !         QuadBound%Xi(1) = - sqrt(3.0_dp / 5.0_dp)
    !         QuadBound%Xi(2) = 0.0_dp
    !         QuadBound%Xi(3) = sqrt(3.0_dp / 5.0_dp)

    !         QuadBound%W(1) = 5.0_dp / 9.0_dp
    !         QuadBound%W(2) = 8.0_dp / 9.0_dp
    !         QuadBound%W(3) = 5.0_dp / 9.0_dp

    !     case(4)
    !         allocate(QuadBound%Xi(4))
    !         allocate(QuadBound%W(4))

    !         QuadBound%Xi(1) = - sqrt((3.0_dp + 2.0_dp * sqrt(6.0_dp / 5.0_dp)) / 7.0_dp)
    !         QuadBound%Xi(2) = - sqrt((3.0_dp - 2.0_dp * sqrt(6.0_dp / 5.0_dp)) / 7.0_dp)
    !         QuadBound%Xi(3) = sqrt((3.0_dp - 2.0_dp * sqrt(6.0_dp / 5.0_dp)) / 7.0_dp)
    !         QuadBound%Xi(4) = sqrt((3.0_dp + 2.0_dp * sqrt(6.0_dp / 5.0_dp)) / 7.0_dp)

    !         QuadBound%W(1) = (18.0_dp - sqrt(30.0_dp)) / 36.0_dp
    !         QuadBound%W(2) = (18.0_dp + sqrt(30.0_dp)) / 36.0_dp
    !         QuadBound%W(3) = (18.0_dp + sqrt(30.0_dp)) / 36.0_dp
    !         QuadBound%W(4) = (18.0_dp - sqrt(30.0_dp)) / 36.0_dp

    !     case(5)
    !         allocate(QuadBound%Xi(5))
    !         allocate(QuadBound%W(5))

    !         QuadBound%Xi(1) = - 0.9061798459386640_dp
    !         QuadBound%Xi(2) = - 0.5384693101056831_dp
    !         QuadBound%Xi(3) = 0.0_dp
    !         QuadBound%Xi(4) = 0.5384693101056831_dp
    !         QuadBound%Xi(5) = 0.9061798459386640_dp

    !         QuadBound%W(1) = 0.2369268850561891_dp
    !         QuadBound%W(2) = 0.4786286704993665_dp
    !         QuadBound%W(3) = 0.5688888888888889_dp
    !         QuadBound%W(4) = 0.4786286704993665_dp
    !         QuadBound%W(5) = 0.2369268850561891_dp
        
    !     case(6)
    !         allocate(QuadBound%Xi(6))
    !         allocate(QuadBound%W(6))

    !         QuadBound%Xi(1) = - 0.9324695142031521_dp
    !         QuadBound%Xi(2) = - 0.6612093864662645_dp
    !         QuadBound%Xi(3) = - 0.2386191860831969_dp
    !         QuadBound%Xi(4) = 0.2386191860831969_dp
    !         QuadBound%Xi(5) = 0.6612093864662645_dp
    !         QuadBound%Xi(6) = 0.9324695142031521_dp

    !         QuadBound%W(1) = 0.1713244923791704_dp
    !         QuadBound%W(2) = 0.3607615730481386_dp
    !         QuadBound%W(3) = 0.4679139345726910_dp
    !         QuadBound%W(4) = 0.4679139345726910_dp
    !         QuadBound%W(5) = 0.3607615730481386_dp
    !         QuadBound%W(6) = 0.1713244923791704_dp

    !     end select
end subroutine TriangleQuadrature

subroutine HexahedralQuadrature(Quad, QuadBound, IntegOrder)

    type(t_Quadrature)  :: Quad, QuadBound
    integer, intent(in) :: IntegOrder
    integer :: ii, j, k, gp
    
    
    allocate(Quad%Xi((IntegOrder+1)**3))
    allocate(Quad%Eta((IntegOrder+1)**3))
    allocate(Quad%Zeta((IntegOrder+1)**3))
    allocate(Quad%W((IntegOrder+1)**3))
    Quad%NoPoints = (IntegOrder+1)**3

    gp = 1
    do ii = 1, IntegOrder + 1
        do j = 1, IntegOrder + 1
            do k = 1, IntegOrder + 1
                Quad%Xi(gp) = QuadBound%Xi(ii)
                Quad%Eta(gp) = QuadBound%Xi(j)
                Quad%Zeta(gp) = QuadBound%Xi(k)
                Quad%W(gp) = QuadBound%W(ii) * QuadBound%W(j) * QuadBound%W(k)
                gp = gp + 1
            end do
        end do
    end do

end subroutine HexahedralQuadrature

subroutine GetQuadTet(Quad)
    type(t_Quadrature)           :: Quad

    select case (Quad%NoPoints)

        case(1)
            allocate(Quad%Xi(1))
            allocate(Quad%Eta(1))
            allocate(Quad%Zeta(1))
            allocate(Quad%W(1))

            Quad%Xi(1) = 1.0_dp / 4.0_dp
            Quad%Eta(1) = 1.0_dp / 4.0_dp
            Quad%Zeta(1) = 1.0_dp / 4.0_dp
            Quad%W(1) = 1.0_dp/6.0_dp

        case(5)
            allocate(Quad%Xi(5))
            allocate(Quad%Eta(5))
            allocate(Quad%Zeta(5))
            allocate(Quad%W(5))

            Quad%Xi(1) = 1.0_dp / 4.0_dp
            Quad%Xi(2) = 1.0_dp / 6.0_dp
            Quad%Xi(3) = 1.0_dp / 6.0_dp
            Quad%Xi(4) = 1.0_dp / 6.0_dp
            Quad%Xi(5) = 1.0_dp / 2.0_dp

            Quad%Eta(1) = 1.0_dp / 4.0_dp
            Quad%Eta(2) = 1.0_dp / 6.0_dp
            Quad%Eta(3) = 1.0_dp / 6.0_dp
            Quad%Eta(4) = 1.0_dp / 2.0_dp
            Quad%Eta(5) = 1.0_dp / 6.0_dp

            Quad%Zeta(1) = 1.0_dp / 4.0_dp
            Quad%Zeta(2) = 1.0_dp / 6.0_dp
            Quad%Zeta(3) = 1.0_dp / 2.0_dp
            Quad%Zeta(4) = 1.0_dp / 6.0_dp
            Quad%Zeta(5) = 1.0_dp / 6.0_dp

            Quad%W(1) = - 2.0_dp / 15.0_dp
            Quad%W(2:5) = 3.0_dp / 40.0_dp

    end select

    ! select case (QuadBound%NoPoints)
    !     case(3)
        
    !         allocate(QuadBound%Xi(3))
    !         allocate(QuadBound%Eta(3))
    !         allocate(QuadBound%Zeta(3))
    !         allocate(QuadBound%W(3))

    !         QuadBound%eta(1) = 0.1666666666666667_8
    !         QuadBound%xi(1) = 0.1666666666666667_8
    !         QuadBound% eta(2) = 0.1666666666666667_8
    !         QuadBound%xi(2) = 0.6666666666666667_8
    !         QuadBound%eta(3) = 0.6666666666666667_8
    !         QuadBound%xi(3) = 0.1666666666666667_8
    !         QuadBound%W(1) = 0.3333333333333333_8
    !         QuadBound%W(2) = 0.3333333333333333_8
    !         QuadBound%W(3) = 0.3333333333333333_8
    !         end select

end subroutine GetQuadTet


subroutine Spectral1DQuadrature(Quad, IntegOrder)
    type(t_Quadrature)      :: Quad
    integer                 :: IntegOrder

    allocate(Quad%Xi(IntegOrder+1))
    allocate(Quad%W(IntegOrder+1))
    Quad%NoPoints = IntegOrder + 1

    Quad%Xi(1) = -1.0_dp
    Quad%Xi(IntegOrder+1) = 1.0_dp

    Quad%W(1) = 2.0_dp / (IntegOrder * (IntegOrder + 1))
    Quad%W(IntegOrder+1) = Quad%W(1)
    
    select case (IntegOrder)
    case(2)
        Quad%Xi(2) = 0.0_dp
        Quad%W(2) = 4.0_dp / 3.0_dp
    case(3)
        Quad%Xi(2) = -1.0_dp / sqrt(5.0_dp)
        Quad%Xi(3) = -Quad%Xi(2)

        Quad%W(2) = 5.0_dp / 6.0_dp
        Quad%W(3) = Quad%W(2)
    case(4)
        Quad%Xi(2) = -sqrt(3.0_dp/7.0_dp)
        Quad%Xi(3) = 0.0_dp
        Quad%Xi(4) = -Quad%Xi(2)

        Quad%W(2) = 49.0_dp / 90.0_dp
        Quad%W(3) = 32.0_dp / 45.0_dp
        Quad%W(4) = Quad%W(2)
    case(5)
        Quad%Xi(2) = -0.765055323929465_dp
        Quad%Xi(3) = -0.285231516480645_dp
        Quad%Xi(4) = -Quad%Xi(3)
        Quad%Xi(5) = -Quad%Xi(2)

        Quad%W(2) = 0.37847495629785_dp
        Quad%W(3) = 0.55485837703549_dp
        Quad%W(4) = Quad%W(3)
        Quad%W(5) = Quad%W(2)
    case(6)
        Quad%Xi(2) = -0.830223896278567_dp
        Quad%Xi(3) = -0.468848793470714_dp
        Quad%Xi(4) = 0.0_dp
        Quad%Xi(5) = -Quad%Xi(3)
        Quad%Xi(6) = -Quad%Xi(2)

        Quad%W(2) = 0.27682604736157_dp
        Quad%W(3) = 0.43174538120986_dp
        Quad%W(4) = 0.48761904761905_dp
        Quad%W(5) = Quad%W(3)
        Quad%W(6) = Quad%W(2)
    case(7)
        Quad%Xi(2) = -0.871740148509606_dp
        Quad%Xi(3) = -0.591700181433142_dp
        Quad%Xi(4) = -0.209299217902478_dp
        Quad%Xi(5) = -Quad%Xi(4)
        Quad%Xi(6) = -Quad%Xi(3)
        Quad%Xi(7) = -Quad%Xi(2)

        Quad%W(2) = 0.21070422714350_dp
        Quad%W(3) = 0.34112269248350_dp
        Quad%W(4) = 0.41245879465870_dp
        Quad%W(5) = Quad%W(4)
        Quad%W(6) = Quad%W(3)
        Quad%W(7) = Quad%W(2)
    case(8)
        Quad%Xi(2) = -0.899757995411460_dp
        Quad%Xi(3) = -0.677186279510737_dp
        Quad%Xi(4) = -0.363117463826178_dp
        Quad%Xi(5) = 0.0_dp
        Quad%Xi(6) = -Quad%Xi(4)
        Quad%Xi(7) = -Quad%Xi(3)
        Quad%Xi(8) = -Quad%Xi(2)

        Quad%W(2) = 0.16549536156080688_dp
        Quad%W(3) = 0.274538712500162_dp
        Quad%W(4) = 0.3464285109730465_dp
        Quad%W(5) = 0.3715192743764172_dp
        Quad%W(6) = Quad%W(4)
        Quad%W(7) = Quad%W(3)
        Quad%W(8) = Quad%W(2)
    case(9)
        Quad%Xi(2) = -0.919533908166459_dp
        Quad%Xi(3) = -0.738773865105505_dp
        Quad%Xi(4) = -0.477924949810444_dp
        Quad%Xi(5) = -0.165278957666387_dp
        Quad%Xi(6) = -Quad%Xi(5)
        Quad%Xi(7) = -Quad%Xi(4)
        Quad%Xi(8) = -Quad%Xi(3)
        Quad%Xi(9) = -Quad%Xi(2)

        Quad%W(2) = 0.13330599085107228_dp
        Quad%W(3) = 0.2248893420631255_dp
        Quad%W(4) = 0.2920426836796838_dp
        Quad%W(5) = 0.32753976118389755_dp
        Quad%W(6) = Quad%W(5)
        Quad%W(7) = Quad%W(4)
        Quad%W(8) = Quad%W(3)
        Quad%W(9) = Quad%W(2) 
    case(10)
        Quad%Xi(2) = -0.934001430408059_dp
        Quad%Xi(3) = -0.784483473663144_dp
        Quad%Xi(4) = -0.565235326996205_dp
        Quad%Xi(5) = -0.295758135586939_dp
        Quad%Xi(6) = 0.0_dp
        Quad%Xi(7) = -Quad%Xi(5)
        Quad%Xi(8) = -Quad%Xi(4)
        Quad%Xi(9) = -Quad%Xi(3)
        Quad%Xi(10) = -Quad%Xi(2)

        Quad%W(2) = 0.10961227326699513_dp
        Quad%W(3) = 0.18716988178030833_dp
        Quad%W(4) = 0.24804810426402857_dp
        Quad%W(5) = 0.2868791247790081_dp
        Quad%W(6) = 0.3002175954556907_dp
        Quad%W(7) = Quad%W(5)
        Quad%W(8) = Quad%W(4)
        Quad%W(9) = Quad%W(3)
        Quad%W(10) = Quad%W(2)

    case(11)
        Quad%Xi(2) = -0.9448992722296681_dp
        Quad%Xi(3) = -0.8192793216440067_dp
        Quad%Xi(4) = -0.6328761530318606_dp
        Quad%Xi(5) = -0.3995309409653489_dp
        Quad%Xi(6) = -0.1365529328549276_dp
        Quad%Xi(7) = -Quad%Xi(6)
        Quad%Xi(8) = -Quad%Xi(5)
        Quad%Xi(9) = -Quad%Xi(4)
        Quad%Xi(10) = -Quad%Xi(3)
        Quad%Xi(11) = -Quad%Xi(2)

        Quad%W(2) = 0.09168451741320352_dp
        Quad%W(3) = 0.15797470556437104_dp
        Quad%W(4) = 0.21250841776102014_dp
        Quad%W(5) = 0.25127560319920128_dp
        Quad%W(6) = 0.2714052409106962_dp
        Quad%W(7) = Quad%W(6)
        Quad%W(8) = Quad%W(5)
        Quad%W(9) = Quad%W(4)
        Quad%W(10) = Quad%W(3)
        Quad%W(11) = Quad%W(2)

    case(12)
        Quad%Xi(2) = -0.9533098466421639_dp
        Quad%Xi(3) = -0.8463475646518723_dp
        Quad%Xi(4) = -0.6861884690817574_dp
        Quad%Xi(5) = -0.4829098210913362_dp
        Quad%Xi(6) = -0.249286930106240_dp
        Quad%Xi(7) = 0.0_dp
        Quad%Xi(8) = -Quad%Xi(6)
        Quad%Xi(9) = -Quad%Xi(5)
        Quad%Xi(10) = -Quad%Xi(4)
        Quad%Xi(11) = -Quad%Xi(3)
        Quad%Xi(12) = -Quad%Xi(2)

        Quad%W(2) = 0.07780168674682487_dp
        Quad%W(3) = 0.13498192668960732_dp
        Quad%W(4) = 0.1836468652035501_dp
        Quad%W(5) = 0.2207677935661101_dp
        Quad%W(6) = 0.2440157903066763_dp
        Quad%W(7) = 0.2519308493334467_dp
        Quad%W(8) = Quad%W(6)
        Quad%W(9) = Quad%W(5)
        Quad%W(10) = Quad%W(4)
        Quad%W(11) = Quad%W(3)
        Quad%W(12) = Quad%W(2)

    case(13)
        Quad%Xi(2) = -0.959935045267261_dp
        Quad%Xi(3) = -0.867801053830347_dp
        Quad%Xi(4) = -0.728868599091326_dp
        Quad%Xi(5) = -0.550639402928647_dp
        Quad%Xi(6) = -0.342724013342712_dp
        Quad%Xi(7) = -0.116331868883703_dp
        Quad%Xi(8) = -Quad%Xi(7)
        Quad%Xi(9) = -Quad%Xi(6)
        Quad%Xi(10) = -Quad%Xi(5)
        Quad%Xi(11) = -Quad%Xi(4)
        Quad%Xi(12) = -Quad%Xi(3)
        Quad%Xi(13) = -Quad%Xi(2)

        Quad%W(2) = 0.06683728449768153_dp
        Quad%W(3) = 0.11658665589871228_dp
        Quad%W(4) = 0.16002185176295067_dp
        Quad%W(5) = 0.1948261493734163_dp
        Quad%W(6) = 0.21912625300977057_dp
        Quad%W(7) = 0.23161279446845698_dp
        Quad%W(8) = Quad%W(7)
        Quad%W(9) = Quad%W(6)
        Quad%W(10) = Quad%W(5)
        Quad%W(11) = Quad%W(4)
        Quad%W(12) = Quad%W(3)
        Quad%W(13) = Quad%W(2)

    case(14)
        Quad%Xi(2) = -0.965245926503839_dp
        Quad%Xi(3) = -0.885082044222976_dp
        Quad%Xi(4) = -0.763519689951815_dp
        Quad%Xi(5) = -0.606253205469845_dp
        Quad%Xi(6) = -0.420638054713672_dp
        Quad%Xi(7) = -0.215353955363794_dp
        Quad%Xi(8) = 0.0_dp
        Quad%Xi(9) = -Quad%Xi(7)
        Quad%Xi(10) = -Quad%Xi(6)
        Quad%Xi(11) = -Quad%Xi(5)
        Quad%Xi(12) = -Quad%Xi(4)
        Quad%Xi(13) = -Quad%Xi(3)
        Quad%Xi(14) = -Quad%Xi(2)

        Quad%W(2) = 0.05802989302860054_dp
        Quad%W(3) = 0.10166007032571654_dp
        Quad%W(4) = 0.1405116998024292_dp
        Quad%W(5) = 0.17278964725360046_dp
        Quad%W(6) = 0.19698723596461332_dp
        Quad%W(7) = 0.211973585926821_dp
        Quad%W(8) = 0.21704811634881566_dp
        Quad%W(9) = Quad%W(7)
        Quad%W(10) = Quad%W(6)
        Quad%W(11) = Quad%W(5)
        Quad%W(12) = Quad%W(4)
        Quad%W(13) = Quad%W(3)
        Quad%W(14) = Quad%W(2)

    case(15)
        Quad%Xi(2) = -0.969568046270218_dp
        Quad%Xi(3) = -0.899200533093472_dp
        Quad%Xi(4) = -0.7920082918618151_dp
        Quad%Xi(5) = -0.652388702882493_dp
        Quad%Xi(6) = -0.486059421887137_dp
        Quad%Xi(7) = -0.2998304689007632_dp
        Quad%Xi(8) = -0.1013262735219491_dp
        Quad%Xi(9) = -Quad%Xi(8)
        Quad%Xi(10) = -Quad%Xi(7)
        Quad%Xi(11) = -Quad%Xi(6)
        Quad%Xi(12) = -Quad%Xi(5)
        Quad%Xi(13) = -Quad%Xi(4)
        Quad%Xi(14) = -Quad%Xi(3)
        Quad%Xi(15) = -Quad%Xi(2)

        Quad%W(2) = 0.05085036100592187_dp
        Quad%W(3) = 0.08939369732593062_dp
        Quad%W(4) = 0.1242553821325135_dp
        Quad%W(5) = 0.15402698080716518_dp
        Quad%W(6) = 0.17749191339170411_dp
        Quad%W(7) = 0.19369002382520362_dp
        Quad%W(8) = 0.20195830817822985_dp
        Quad%W(9) = 0.20195830817822985_dp
        Quad%W(10) = Quad%W(7)
        Quad%W(11) = Quad%W(6)
        Quad%W(12) = Quad%W(5)
        Quad%W(13) = Quad%W(4)
        Quad%W(14) = Quad%W(3)
        Quad%W(15) = Quad%W(2)

    case(16)
        Quad%Xi(2) = -0.9731321766314183_dp
        Quad%Xi(3) = -0.910879995915574_dp
        Quad%Xi(4) = -0.8156962512217703_dp
        Quad%Xi(5) = -0.6910289806276847_dp
        Quad%Xi(6) = -0.541385399330102_dp
        Quad%Xi(7) = -0.3721744335654772_dp
        Quad%Xi(8) = -0.189511973518317_dp
        Quad%Xi(9) = 0.0_dp
        Quad%Xi(10) = -Quad%Xi(8)
        Quad%Xi(11) = -Quad%Xi(7)
        Quad%Xi(12) = -Quad%Xi(6)
        Quad%Xi(13) = -Quad%Xi(5)
        Quad%Xi(14) = -Quad%Xi(4)
        Quad%Xi(15) = -Quad%Xi(3)
        Quad%Xi(16) = -Quad%Xi(2)

        Quad%W(2) = 0.04492194054325292_dp
        Quad%W(3) = 0.07919827050368623_dp
        Quad%W(4) = 0.11059290900702798_dp
        Quad%W(5) = 0.13798774620192722_dp
        Quad%W(6) = 0.1603946619976215_dp
        Quad%W(7) = 0.1770042535156577_dp
        Quad%W(8) = 0.18721633967761928_dp
        Quad%W(9) = 0.19066187475346943_dp
        Quad%W(10) = Quad%W(8)
        Quad%W(11) = Quad%W(7)
        Quad%W(12) = Quad%W(6)
        Quad%W(13) = Quad%W(5)
        Quad%W(14) = Quad%W(4)
        Quad%W(15) = Quad%W(3)
        Quad%W(16) = Quad%W(2)

    case(17)
        Quad%Xi(2) = -0.976105557412198_dp
        Quad%Xi(3) = -0.920649185347533_dp
        Quad%Xi(4) = -0.835593535218090_dp
        Quad%Xi(5) = -0.723679329283243_dp
        Quad%Xi(6) = -0.588504834318661_dp
        Quad%Xi(7) = -0.434415036912123_dp
        Quad%Xi(8) = -0.2663626528782805_dp
        Quad%Xi(9) = -0.089749093484652_dp
        Quad%Xi(10) = -Quad%Xi(9)
        Quad%Xi(11) = -Quad%Xi(8)
        Quad%Xi(12) = -Quad%Xi(7)
        Quad%Xi(13) = -Quad%Xi(6)
        Quad%Xi(14) = -Quad%Xi(5)
        Quad%Xi(15) = -Quad%Xi(4)
        Quad%Xi(16) = -Quad%Xi(3)
        Quad%Xi(17) = -Quad%Xi(2)

        Quad%W(2) = 0.039970628810914184_dp
        Quad%W(3) = 0.07063716688563393_dp
        Quad%W(4) = 0.0990162717175025_dp
        Quad%W(5) = 0.12421053313296582_dp
        Quad%W(6) = 0.1454119615738022_dp
        Quad%W(7) = 0.16193951723760272_dp
        Quad%W(8) = 0.17326210948945625_dp
        Quad%W(9) = 0.17901586343970305_dp
        Quad%W(10) = 0.17901586343970305_dp
        Quad%W(11) = Quad%W(8)
        Quad%W(12) = Quad%W(7)
        Quad%W(13) = Quad%W(6)
        Quad%W(14) = Quad%W(5)
        Quad%W(15) = Quad%W(4)
        Quad%W(16) = Quad%W(3)
        Quad%W(17) = Quad%W(2)

    case(18)
        Quad%Xi(2) = -0.978611766222080_dp
        Quad%Xi(3) = -0.928901528152586_dp
        Quad%Xi(4) = -0.852460577796646_dp
        Quad%Xi(5) = -0.751494202552613_dp
        Quad%Xi(6) = -0.628908137265221_dp
        Quad%Xi(7) = -0.488229285680714_dp
        Quad%Xi(8) = -0.333504847824499_dp
        Quad%Xi(9) = -0.169186023409282_dp
        Quad%Xi(10) = 0.0_dp
        Quad%Xi(11) = -Quad%Xi(9)
        Quad%Xi(12) = -Quad%Xi(8)
        Quad%Xi(13) = -Quad%Xi(7)
        Quad%Xi(14) = -Quad%Xi(6)
        Quad%Xi(15) = -Quad%Xi(5)
        Quad%Xi(16) = -Quad%Xi(4)
        Quad%Xi(17) = -Quad%Xi(3)
        Quad%Xi(18) = -Quad%Xi(2)

        Quad%W(2) = 0.035793365186175874_dp
        Quad%W(3) = 0.0633818917626272_dp
        Quad%W(4) = 0.08913175709920798_dp
        Quad%W(5) = 0.11231534147730572_dp
        Quad%W(6) = 0.1322672804487499_dp
        Quad%W(7) = 0.14841394259593893_dp
        Quad%W(8) = 0.16029092404406128_dp
        Quad%W(9) = 0.16755658452714284_dp
        Quad%W(10) = 0.17000191928482725_dp
        Quad%W(11) = Quad%W(9)
        Quad%W(12) = Quad%W(8)
        Quad%W(13) = Quad%W(7)
        Quad%W(14) = Quad%W(6)
        Quad%W(15) = Quad%W(5)
        Quad%W(16) = Quad%W(4)
        Quad%W(17) = Quad%W(3)
        Quad%W(18) = Quad%W(2)

    case(19)
        Quad%Xi(2) = -0.980743704893914_dp
        Quad%Xi(3) = -0.935934498812665_dp
        Quad%Xi(4) = -0.866877978089950_dp
        Quad%Xi(5) = -0.775368260952056_dp
        Quad%Xi(6) = -0.663776402290311_dp
        Quad%Xi(7) = -0.534992864031886_dp
        Quad%Xi(8) = -0.392353183713909_dp
        Quad%Xi(9) = -0.239551705922986_dp
        Quad%Xi(10) = -0.080545937238822_dp
        Quad%Xi(11) = -Quad%Xi(10)
        Quad%Xi(12) = -Quad%Xi(9)
        Quad%Xi(13) = -Quad%Xi(8)
        Quad%Xi(14) = -Quad%Xi(7)
        Quad%Xi(15) = -Quad%Xi(6)
        Quad%Xi(16) = -Quad%Xi(5)
        Quad%Xi(17) = -Quad%Xi(4)
        Quad%Xi(18) = -Quad%Xi(3)
        Quad%Xi(19) = -Quad%Xi(2)

        Quad%W(2) = 0.03223712318848816_dp
        Quad%W(3) = 0.05718180212756649_dp
        Quad%W(4) = 0.08063176399612001_dp
        Quad%W(5) = 0.10199149969945108_dp
        Quad%W(6) = 0.12070922762867593_dp
        Quad%W(7) = 0.1363004823587244_dp
        Quad%W(8) = 0.1483615540709169_dp
        Quad%W(9) = 0.15658010264747546_dp
        Quad%W(10) = 0.16074328638784577_dp
        Quad%W(11) = 0.16074328638784577_dp
        Quad%W(12) = Quad%W(9)
        Quad%W(13) = Quad%W(8)
        Quad%W(14) = Quad%W(7)
        Quad%W(15) = Quad%W(6)
        Quad%W(16) = Quad%W(5)
        Quad%W(17) = Quad%W(4)
        Quad%W(18) = Quad%W(3)
        Quad%W(19) = Quad%W(2)

    case(20)
        Quad%Xi(2) = -0.982572296604548_dp
        Quad%Xi(3) = -0.941976296959746_dp
        Quad%Xi(4) = -0.879294755323591_dp
        Quad%Xi(5) = -0.796001926077712_dp
        Quad%Xi(6) = -0.694051026062223_dp
        Quad%Xi(7) = -0.575831960261831_dp
        Quad%Xi(8) = -0.444115783279002_dp
        Quad%Xi(9) = -0.301989856508765_dp
        Quad%Xi(10) = -0.152785515802186_dp
        Quad%Xi(11) = 0.0_dp
        Quad%Xi(12) = -Quad%Xi(10)
        Quad%Xi(13) = -Quad%Xi(9)
        Quad%Xi(14) = -Quad%Xi(8)
        Quad%Xi(15) = -Quad%Xi(7)
        Quad%Xi(16) = -Quad%Xi(6)
        Quad%Xi(17) = -Quad%Xi(5)
        Quad%Xi(18) = -Quad%Xi(4)
        Quad%Xi(19) = -Quad%Xi(3)
        Quad%Xi(20) = -Quad%Xi(2)

        Quad%W(2) = 0.029184840098506866_dp
        Quad%W(3) = 0.05184316900084789_dp
        Quad%W(4) = 0.07327391818507369_dp
        Quad%W(5) = 0.09298546795788497_dp
        Quad%W(6) = 0.1105170832191237_dp
        Quad%W(7) = 0.12545812119086924_dp
        Quad%W(8) = 0.13745846286004137_dp
        Quad%W(9) = 0.14623686244797748_dp
        Quad%W(10) = 0.1515875751116814_dp
        Quad%W(11) = 0.15338519033217496_dp
        Quad%W(12) = Quad%W(10)
        Quad%W(13) = Quad%W(9)
        Quad%W(14) = Quad%W(8)
        Quad%W(15) = Quad%W(7)
        Quad%W(16) = Quad%W(6)
        Quad%W(17) = Quad%W(5)
        Quad%W(18) = Quad%W(4)
        Quad%W(19) = Quad%W(3)
        Quad%W(20) = Quad%W(2)
        
    end select

end subroutine Spectral1DQuadrature

subroutine Spectral2DQuadrature(Quad, QuadBound)
    type(t_Quadrature)  :: Quad
    type(t_Quadrature)  :: QuadBound
    integer                             :: ii, j, k

    Quad%NoPoints = QuadBound%NoPoints ** 2

    allocate(Quad%Xi(Quad%NoPoints))
    allocate(Quad%Eta(Quad%NoPoints))
    allocate(Quad%W(Quad%NoPoints))

    k = 1
    do ii = 1, QuadBound%NoPoints
        do j = 1, QuadBound%NoPoints
            Quad%Xi(k) = QuadBound%Xi(ii)
            Quad%Eta(k) = QuadBound%Xi(j)
            Quad%W(k) = QuadBound%W(ii) * QuadBound%W(j)
            k = k + 1
        end do
    end do
    
end subroutine Spectral2DQuadrature

end module m_quadrature
