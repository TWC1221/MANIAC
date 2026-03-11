module m_quadrature
!-----------------------------------------------------------------------!
! Purpose:                                                            -!
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
! Record of revisions:                                                -!
!   Date       Programmer     Description of change                    -!
!   ====       ==========     =====================                    -!
! 14/03/24      C. Jones         Original code                         -!    
! 12/02/2026    T. Charlton      Quad Generation via LaPack             -!
!-----------------------------------------------------------------------!

use m_constants
implicit none

type :: t_quadrature
    real(dp), dimension(:), allocatable :: Xi, Eta, Zeta, W
    real(dp), allocatable :: dir(:,:)
    integer :: NoPoints
end type t_quadrature

contains

subroutine Get1DLineQuad(Quad1D, IntegOrder, flag_adjoint)
    type(t_quadrature) :: Quad1D
    integer, intent(in) :: IntegOrder
    logical, intent(in), optional :: flag_adjoint
    
    real(dp), allocatable :: d(:), e(:), z(:,:)
    integer :: info
    real(dp), allocatable :: work(:)
    integer :: i

    Quad1D%NoPoints = IntegOrder
    if (allocated(Quad1D%Xi)) deallocate(Quad1D%Xi)
    if (allocated(Quad1D%W))  deallocate(Quad1D%W)
    allocate(Quad1D%Xi(IntegOrder), Quad1D%W(IntegOrder))
    allocate(d(IntegOrder), e(IntegOrder-1), z(IntegOrder, IntegOrder))
    
    d = 0.0_dp
    do i = 1, IntegOrder - 1
        e(i) = real(i, dp) / sqrt(4.0_dp * i**2 - 1.0_dp)
    end do

    z = 0.0_dp
    do i = 1, IntegOrder
        z(i,i) = 1.0_dp
    end do

    allocate(work(max(1, 2*IntegOrder-2)))

    call DSTEQR('V', IntegOrder, d, e, z, IntegOrder, work, info)

    if (info /= 0) stop "Error: Eigensolver failed in Get1DLineQuad"

    Quad1D%Xi = d
    do i = 1, IntegOrder
        Quad1D%W(i) = 2.0_dp * z(1, i)**2
    end do

    if (present(flag_adjoint)) Quad1D%Xi = - Quad1D%Xi

end subroutine Get1DLineQuad

subroutine HexahedralQuadrature(Quad, IntegOrder, writeout)
    type(t_Quadrature), intent(out)  :: Quad
    integer, intent(in) :: IntegOrder
    logical, intent(in) :: writeout
    type(t_Quadrature)  :: Quad1D
    integer :: i, j, k, gp
    
    call Get1DLineQuad(Quad1D, IntegOrder)

    if (writeout) write(*,*) ">>> Hexahedral Quadrature Master Element Points ..."
    Quad%NoPoints = IntegOrder**3
    allocate(Quad%Xi(Quad%NoPoints), Quad%Eta(Quad%NoPoints), Quad%Zeta(Quad%NoPoints), Quad%W(Quad%NoPoints))

    gp = 1
    do k = 1, IntegOrder
        do j = 1, IntegOrder
            do i = 1, IntegOrder
                Quad%Xi(gp)   = Quad1D%Xi(i)
                Quad%Eta(gp)  = Quad1D%Xi(j)
                Quad%Zeta(gp) = Quad1D%Xi(k)
                Quad%W(gp)    = Quad1D%W(i) * Quad1D%W(j) * Quad1D%W(k)

                if (writeout) write(*,*) Quad%Xi(gp), Quad%Eta(gp), Quad%Zeta(gp), Quad%W(gp)
                gp = gp + 1
            end do
        end do
    end do

end subroutine HexahedralQuadrature

subroutine QuadrilateralQuadrature(Quad, IntegOrder)
    type(t_quadrature), intent(out) :: Quad
    integer, intent(in) :: IntegOrder
    type(t_quadrature) :: Quad1D
    integer :: i, j, counter

    call Get1DLineQuad(Quad1D, IntegOrder)

    Quad%NoPoints = IntegOrder**2
    allocate(Quad%Xi(Quad%NoPoints), Quad%Eta(Quad%NoPoints), Quad%W(Quad%NoPoints))
    if (allocated(Quad%Zeta)) deallocate(Quad%Zeta)

    counter = 1
    do j = 1, IntegOrder
        do i = 1, IntegOrder
            Quad%Xi(counter)  = Quad1D%Xi(i)
            Quad%Eta(counter) = Quad1D%Xi(j)
            Quad%W(counter)   = Quad1D%W(i) * Quad1D%W(j)
            counter = counter + 1
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
    integer                             :: ii, jj, k

    Quad%NoPoints = QuadBound%NoPoints ** 2

    allocate(Quad%Xi(Quad%NoPoints))
    allocate(Quad%Eta(Quad%NoPoints))
    allocate(Quad%W(Quad%NoPoints))

    k = 1
    do ii = 1, QuadBound%NoPoints
        do jj = 1, QuadBound%NoPoints
            Quad%Xi(k) = QuadBound%Xi(ii)
            Quad%Eta(k) = QuadBound%Xi(jj)
            Quad%W(k) = QuadBound%W(ii) * QuadBound%W(jj)
            k = k + 1
        end do
    end do
    
end subroutine Spectral2DQuadrature

end module m_quadrature
