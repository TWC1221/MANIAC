module m_quadrature
!-----------------------------------------------------------------------!
! Purpose:                                                             -!
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
! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                    -!
!   ====       ==========     =====================                    -!
! 14/03/24      C. Jones         Original code                         -!    
! 12/02/2026    T. Charlton      Quad Generation via LAPACK            -!
! 16/03/2026    T. Charlton      Angular Generation Algo added         -!
!-----------------------------------------------------------------------!

use m_constants
use m_types
implicit none

public :: t_quadrature, t_sn_quadrature

contains

subroutine InitialiseQuadrature(FE, mesh, Quad1D, Quad, QuadSn, is_adjoint, boundary_order_plus)
    type(t_quadrature), intent(inout)       :: Quad1D, Quad
    type(t_sn_quadrature), intent(inout)    :: QuadSn
    type(t_finite), intent(inout)              :: FE   
    type(t_mesh), intent(in)                :: mesh
    logical, intent(in)                     :: is_adjoint
    integer, intent(in), optional           :: boundary_order_plus
    integer :: b_order, v_order
    type(t_quadrature) :: Quad1D_vol

    v_order = FE%order + 1
    b_order = v_order
    if (present(boundary_order_plus)) b_order = b_order + boundary_order_plus

    call LinearQuadrature(Quad1D, b_order)
    call LinearQuadrature(Quad1D_vol, v_order)
    call QuadrilateralQuadrature(Quad, v_order)

    call AngularQuadrature(mesh%dim, QuadSn%order, QuadSn, is_adjoint)

end subroutine InitialiseQuadrature

subroutine LinearQuadrature(Quad1D, IntegOrder)
    type(t_quadrature) :: Quad1D
    integer, intent(in) :: IntegOrder
    
    real(dp), allocatable :: d(:), e(:), z(:,:)
    integer :: info
    real(dp), allocatable :: work(:)
    integer :: i

    Quad1D%n_points = IntegOrder
    if (allocated(Quad1D%xi)) deallocate(Quad1D%xi)
    if (allocated(Quad1D%weights))  deallocate(Quad1D%weights)
    allocate(Quad1D%xi(IntegOrder), Quad1D%weights(IntegOrder))
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

    Quad1D%xi = d
    do i = 1, IntegOrder
        Quad1D%weights(i) = 2.0_dp * z(1, i)**2
    end do

end subroutine LinearQuadrature

subroutine QuadrilateralQuadrature(Quad, IntegOrder)
    type(t_quadrature), intent(out) :: Quad
    integer, intent(in) :: IntegOrder
    type(t_quadrature) :: Quad1D
    integer :: i, j, counter = 1

    call LinearQuadrature(Quad1D, IntegOrder)

    Quad%n_points = IntegOrder**2
    allocate(Quad%xi(Quad%n_points), Quad%eta(Quad%n_points), Quad%weights(Quad%n_points))

    do i = 1, IntegOrder
        do j = 1, IntegOrder
            Quad%xi(counter)  = Quad1D%xi(i)
            Quad%eta(counter) = Quad1D%xi(j)
            Quad%weights(counter)   = Quad1D%weights(i) * Quad1D%weights(j)
            counter = counter + 1
        end do
    end do
    
end subroutine QuadrilateralQuadrature

subroutine AngularQuadrature(dim, SN, QuadSn, Adjoint)
    integer, intent(in)            :: dim, SN
    type(t_sn_quadrature), intent(inout) :: QuadSn
    logical, intent(in)            :: Adjoint
    
    integer  :: i, j, k, m, n_levels, n_octant, q, i_sign
    integer  :: ids(3)
    real(dp), allocatable :: levels(:)
    real(dp) :: signs(8, 3), m0, m2
    
    n_levels = SN / 2
    n_octant = (n_levels * (n_levels + 1)) / 2
    QuadSn%n_angles = n_octant * (merge(4, 8, dim == 2))

    if (allocated(QuadSn%dirs)) deallocate(QuadSn%dirs)
    if (allocated(QuadSn%weights))      deallocate(QuadSn%weights)
    allocate(QuadSn%dirs(QuadSn%n_angles, 3), QuadSn%weights(QuadSn%n_angles))
    allocate(levels(n_levels))

    select case(SN)
        case(2);  levels = [0.57735027_dp]
        case(4);  levels = [0.35002120_dp, 0.86889030_dp]
        case(6);  levels = [0.26663550_dp, 0.68150760_dp, 0.92618080_dp]
        case(8);  levels = [0.21821790_dp, 0.57735030_dp, 0.78679580_dp, 0.95118970_dp]
        case(12); levels = [0.16721260_dp, 0.45954760_dp, 0.62802360_dp, 0.76002100_dp, 0.87227060_dp, 0.97163770_dp]
        case(16); levels = [0.13895680_dp, 0.39228930_dp, 0.53709660_dp, 0.65042640_dp, 0.74675060_dp, 0.83199660_dp, 0.90928550_dp, 0.98050090_dp]
        case default; stop "AngularQuadrature: SN order not supported."
    end select

    m = 0
    do i = 1, n_levels
        do j = 1, n_levels - i + 1
            k = n_levels - i - j + 2
            m = m + 1
            QuadSn%dirs(m, :) = [levels(i), levels(j), levels(k)]
            
            ids = [i, j, k]
            if (ids(1) > ids(2)) call swap(ids(1), ids(2))
            if (ids(2) > ids(3)) call swap(ids(2), ids(3))
            if (ids(1) > ids(2)) call swap(ids(1), ids(2))

            select case(SN)
                case(2, 4); QuadSn%weights(m) = 1.0_dp / 3.0_dp
                case(6)
                    if (ids(1)==1 .and. ids(2)==1) then; QuadSn%weights(m) = 0.1761263_dp
                    else;                                QuadSn%weights(m) = 0.1572071_dp; endif
                case(8)
                    if (ids(1)==1 .and. ids(2)==1) then;      QuadSn%weights(m) = 0.1209877_dp
                    else if (ids(1)==2 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0925926_dp
                    else;                                     QuadSn%weights(m) = 0.0907407_dp; endif
                case(12)
                    if (ids(1)==1 .and. ids(2)==1) then;      QuadSn%weights(m) = 0.0707626_dp
                    else if (ids(1)==1 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0558811_dp
                    else if (ids(1)==1 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0373377_dp
                    else if (ids(1)==2 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0502819_dp
                    else;                                     QuadSn%weights(m) = 0.0258513_dp; endif
                case(16)
                    if (ids(1)==1 .and. ids(2)==1) then;      QuadSn%weights(m) = 0.0489872_dp
                    else if (ids(1)==1 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0413296_dp
                    else if (ids(1)==1 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0212326_dp
                    else if (ids(1)==1 .and. ids(2)==4) then; QuadSn%weights(m) = 0.0256207_dp
                    else if (ids(1)==2 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0360486_dp
                    else if (ids(1)==2 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0144589_dp
                    else if (ids(1)==2 .and. ids(2)==4) then; QuadSn%weights(m) = 0.0344958_dp
                    else if (ids(1)==3 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0085179_dp
                    else if (ids(1)==3 .and. ids(2)==4) then; QuadSn%weights(m) = 0.0144589_dp
                    else;                                     QuadSn%weights(m) = 0.0256207_dp; endif
            end select
        end do
    end do

    signs(1,:) = [ 1.0,  1.0,  1.0]; signs(2,:) = [-1.0,  1.0,  1.0]
    signs(3,:) = [-1.0, -1.0,  1.0]; signs(4,:) = [ 1.0, -1.0,  1.0]
    signs(5,:) = [ 1.0,  1.0, -1.0]; signs(6,:) = [-1.0,  1.0, -1.0]
    signs(7,:) = [-1.0, -1.0, -1.0]; signs(8,:) = [ 1.0, -1.0, -1.0]

    do q = 2, (merge(4, 8, dim == 2))
        do i_sign = 1, n_octant
            m = m + 1
            QuadSn%dirs(m, :) = QuadSn%dirs(i_sign, :) * signs(q, :)
            QuadSn%weights(m) = QuadSn%weights(i_sign)
        end do
    end do

    ! Apply Adjoint reversal to ALL angles after generation
    if (Adjoint) then
        QuadSn%dirs = -QuadSn%dirs
    end if

    QuadSn%weights = QuadSn%weights / sum(QuadSn%weights)
    
    m0 = sum(QuadSn%weights)
    m2 = sum(QuadSn%weights * QuadSn%dirs(:,1)**2)
    write(*,'(A,I2,A,F10.8,A,F10.8)') "S", SN, " Verification: M0=", m0, " M2=", m2*3.0_dp

    contains
        subroutine swap(a, b)
            integer, intent(inout) :: a, b
            integer :: tmp
            tmp = a; a = b; b = tmp
        end subroutine
end subroutine AngularQuadrature

subroutine HexahedralQuadrature(Quad, QuadBound, IntegOrder)

    type(t_Quadrature)  :: Quad, QuadBound
    integer, intent(in) :: IntegOrder
    integer :: ii, jj, k, gp
    
    
    allocate(Quad%xi((IntegOrder+1)**3))
    allocate(Quad%eta((IntegOrder+1)**3))
    allocate(Quad%zeta((IntegOrder+1)**3))
    allocate(Quad%weights((IntegOrder+1)**3))
    Quad%n_points = (IntegOrder+1)**3

    gp = 1
    do ii = 1, IntegOrder + 1
        do jj = 1, IntegOrder + 1
            do k = 1, IntegOrder + 1
                Quad%xi(gp) = QuadBound%xi(ii)
                Quad%eta(gp) = QuadBound%xi(jj)
                Quad%zeta(gp) = QuadBound%xi(k)
                Quad%weights(gp) = QuadBound%weights(ii) * QuadBound%weights(jj) * QuadBound%weights(k)
                gp = gp + 1
            end do
        end do
    end do

end subroutine HexahedralQuadrature

subroutine Spectral1DQuadrature(Quad, IntegOrder)
    type(t_Quadrature)      :: Quad
    integer                 :: IntegOrder

    allocate(Quad%xi(IntegOrder+1))
    allocate(Quad%weights(IntegOrder+1))
    Quad%n_points = IntegOrder + 1

    Quad%xi(1) = -1.0_dp
    Quad%xi(IntegOrder+1) = 1.0_dp

    Quad%weights(1) = 2.0_dp / (IntegOrder * (IntegOrder + 1))
    Quad%weights(IntegOrder+1) = Quad%weights(1)
    
    select case (IntegOrder)
    case(2)
        Quad%xi(2) = 0.0_dp
        Quad%weights(2) = 4.0_dp / 3.0_dp
    case(3)
        Quad%xi(2) = -1.0_dp / sqrt(5.0_dp)
        Quad%xi(3) = -Quad%xi(2)

        Quad%weights(2) = 5.0_dp / 6.0_dp
        Quad%weights(3) = Quad%weights(2)
    case(4)
        Quad%xi(2) = -sqrt(3.0_dp/7.0_dp)
        Quad%xi(3) = 0.0_dp
        Quad%xi(4) = -Quad%xi(2)

        Quad%weights(2) = 49.0_dp / 90.0_dp
        Quad%weights(3) = 32.0_dp / 45.0_dp
        Quad%weights(4) = Quad%weights(2)
    case(5)
        Quad%xi(2) = -0.765055323929465_dp
        Quad%xi(3) = -0.285231516480645_dp
        Quad%xi(4) = -Quad%xi(3)
        Quad%xi(5) = -Quad%xi(2)

        Quad%weights(2) = 0.37847495629785_dp
        Quad%weights(3) = 0.55485837703549_dp
        Quad%weights(4) = Quad%weights(3)
        Quad%weights(5) = Quad%weights(2)
    case(6)
        Quad%xi(2) = -0.830223896278567_dp
        Quad%xi(3) = -0.468848793470714_dp
        Quad%xi(4) = 0.0_dp
        Quad%xi(5) = -Quad%xi(3)
        Quad%xi(6) = -Quad%xi(2)

        Quad%weights(2) = 0.27682604736157_dp
        Quad%weights(3) = 0.43174538120986_dp
        Quad%weights(4) = 0.48761904761905_dp
        Quad%weights(5) = Quad%weights(3)
        Quad%weights(6) = Quad%weights(2)
    case(7)
        Quad%xi(2) = -0.871740148509606_dp
        Quad%xi(3) = -0.591700181433142_dp
        Quad%xi(4) = -0.209299217902478_dp
        Quad%xi(5) = -Quad%xi(4)
        Quad%xi(6) = -Quad%xi(3)
        Quad%xi(7) = -Quad%xi(2)

        Quad%weights(2) = 0.21070422714350_dp
        Quad%weights(3) = 0.34112269248350_dp
        Quad%weights(4) = 0.41245879465870_dp
        Quad%weights(5) = Quad%weights(4)
        Quad%weights(6) = Quad%weights(3)
        Quad%weights(7) = Quad%weights(2)
    case(8)
        Quad%xi(2) = -0.899757995411460_dp
        Quad%xi(3) = -0.677186279510737_dp
        Quad%xi(4) = -0.363117463826178_dp
        Quad%xi(5) = 0.0_dp
        Quad%xi(6) = -Quad%xi(4)
        Quad%xi(7) = -Quad%xi(3)
        Quad%xi(8) = -Quad%xi(2)

        Quad%weights(2) = 0.16549536156080688_dp
        Quad%weights(3) = 0.274538712500162_dp
        Quad%weights(4) = 0.3464285109730465_dp
        Quad%weights(5) = 0.3715192743764172_dp
        Quad%weights(6) = Quad%weights(4)
        Quad%weights(7) = Quad%weights(3)
        Quad%weights(8) = Quad%weights(2)
    case(9)
        Quad%xi(2) = -0.919533908166459_dp
        Quad%xi(3) = -0.738773865105505_dp
        Quad%xi(4) = -0.477924949810444_dp
        Quad%xi(5) = -0.165278957666387_dp
        Quad%xi(6) = -Quad%xi(5)
        Quad%xi(7) = -Quad%xi(4)
        Quad%xi(8) = -Quad%xi(3)
        Quad%xi(9) = -Quad%xi(2)

        Quad%weights(2) = 0.13330599085107228_dp
        Quad%weights(3) = 0.2248893420631255_dp
        Quad%weights(4) = 0.2920426836796838_dp
        Quad%weights(5) = 0.32753976118389755_dp
        Quad%weights(6) = Quad%weights(5)
        Quad%weights(7) = Quad%weights(4)
        Quad%weights(8) = Quad%weights(3)
        Quad%weights(9) = Quad%weights(2) 
    case(10)
        Quad%xi(2) = -0.934001430408059_dp
        Quad%xi(3) = -0.784483473663144_dp
        Quad%xi(4) = -0.565235326996205_dp
        Quad%xi(5) = -0.295758135586939_dp
        Quad%xi(6) = 0.0_dp
        Quad%xi(7) = -Quad%xi(5)
        Quad%xi(8) = -Quad%xi(4)
        Quad%xi(9) = -Quad%xi(3)
        Quad%xi(10) = -Quad%xi(2)

        Quad%weights(2) = 0.10961227326699513_dp
        Quad%weights(3) = 0.18716988178030833_dp
        Quad%weights(4) = 0.24804810426402857_dp
        Quad%weights(5) = 0.2868791247790081_dp
        Quad%weights(6) = 0.3002175954556907_dp
        Quad%weights(7) = Quad%weights(5)
        Quad%weights(8) = Quad%weights(4)
        Quad%weights(9) = Quad%weights(3)
        Quad%weights(10) = Quad%weights(2)

    case(11)
        Quad%xi(2) = -0.9448992722296681_dp
        Quad%xi(3) = -0.8192793216440067_dp
        Quad%xi(4) = -0.6328761530318606_dp
        Quad%xi(5) = -0.3995309409653489_dp
        Quad%xi(6) = -0.1365529328549276_dp
        Quad%xi(7) = -Quad%xi(6)
        Quad%xi(8) = -Quad%xi(5)
        Quad%xi(9) = -Quad%xi(4)
        Quad%xi(10) = -Quad%xi(3)
        Quad%xi(11) = -Quad%xi(2)

        Quad%weights(2) = 0.09168451741320352_dp
        Quad%weights(3) = 0.15797470556437104_dp
        Quad%weights(4) = 0.21250841776102014_dp
        Quad%weights(5) = 0.25127560319920128_dp
        Quad%weights(6) = 0.2714052409106962_dp
        Quad%weights(7) = Quad%weights(6)
        Quad%weights(8) = Quad%weights(5)
        Quad%weights(9) = Quad%weights(4)
        Quad%weights(10) = Quad%weights(3)
        Quad%weights(11) = Quad%weights(2)

    case(12)
        Quad%xi(2) = -0.9533098466421639_dp
        Quad%xi(3) = -0.8463475646518723_dp
        Quad%xi(4) = -0.6861884690817574_dp
        Quad%xi(5) = -0.4829098210913362_dp
        Quad%xi(6) = -0.249286930106240_dp
        Quad%xi(7) = 0.0_dp
        Quad%xi(8) = -Quad%xi(6)
        Quad%xi(9) = -Quad%xi(5)
        Quad%xi(10) = -Quad%xi(4)
        Quad%xi(11) = -Quad%xi(3)
        Quad%xi(12) = -Quad%xi(2)

        Quad%weights(2) = 0.07780168674682487_dp
        Quad%weights(3) = 0.13498192668960732_dp
        Quad%weights(4) = 0.1836468652035501_dp
        Quad%weights(5) = 0.2207677935661101_dp
        Quad%weights(6) = 0.2440157903066763_dp
        Quad%weights(7) = 0.2519308493334467_dp
        Quad%weights(8) = Quad%weights(6)
        Quad%weights(9) = Quad%weights(5)
        Quad%weights(10) = Quad%weights(4)
        Quad%weights(11) = Quad%weights(3)
        Quad%weights(12) = Quad%weights(2)

    case(13)
        Quad%xi(2) = -0.959935045267261_dp
        Quad%xi(3) = -0.867801053830347_dp
        Quad%xi(4) = -0.728868599091326_dp
        Quad%xi(5) = -0.550639402928647_dp
        Quad%xi(6) = -0.342724013342712_dp
        Quad%xi(7) = -0.116331868883703_dp
        Quad%xi(8) = -Quad%xi(7)
        Quad%xi(9) = -Quad%xi(6)
        Quad%xi(10) = -Quad%xi(5)
        Quad%xi(11) = -Quad%xi(4)
        Quad%xi(12) = -Quad%xi(3)
        Quad%xi(13) = -Quad%xi(2)

        Quad%weights(2) = 0.06683728449768153_dp
        Quad%weights(3) = 0.11658665589871228_dp
        Quad%weights(4) = 0.16002185176295067_dp
        Quad%weights(5) = 0.1948261493734163_dp
        Quad%weights(6) = 0.21912625300977057_dp
        Quad%weights(7) = 0.23161279446845698_dp
        Quad%weights(8) = Quad%weights(7)
        Quad%weights(9) = Quad%weights(6)
        Quad%weights(10) = Quad%weights(5)
        Quad%weights(11) = Quad%weights(4)
        Quad%weights(12) = Quad%weights(3)
        Quad%weights(13) = Quad%weights(2)

    case(14)
        Quad%xi(2) = -0.965245926503839_dp
        Quad%xi(3) = -0.885082044222976_dp
        Quad%xi(4) = -0.763519689951815_dp
        Quad%xi(5) = -0.606253205469845_dp
        Quad%xi(6) = -0.420638054713672_dp
        Quad%xi(7) = -0.215353955363794_dp
        Quad%xi(8) = 0.0_dp
        Quad%xi(9) = -Quad%xi(7)
        Quad%xi(10) = -Quad%xi(6)
        Quad%xi(11) = -Quad%xi(5)
        Quad%xi(12) = -Quad%xi(4)
        Quad%xi(13) = -Quad%xi(3)
        Quad%xi(14) = -Quad%xi(2)

        Quad%weights(2) = 0.05802989302860054_dp
        Quad%weights(3) = 0.10166007032571654_dp
        Quad%weights(4) = 0.1405116998024292_dp
        Quad%weights(5) = 0.17278964725360046_dp
        Quad%weights(6) = 0.19698723596461332_dp
        Quad%weights(7) = 0.211973585926821_dp
        Quad%weights(8) = 0.21704811634881566_dp
        Quad%weights(9) = Quad%weights(7)
        Quad%weights(10) = Quad%weights(6)
        Quad%weights(11) = Quad%weights(5)
        Quad%weights(12) = Quad%weights(4)
        Quad%weights(13) = Quad%weights(3)
        Quad%weights(14) = Quad%weights(2)

    case(15)
        Quad%xi(2) = -0.969568046270218_dp
        Quad%xi(3) = -0.899200533093472_dp
        Quad%xi(4) = -0.7920082918618151_dp
        Quad%xi(5) = -0.652388702882493_dp
        Quad%xi(6) = -0.486059421887137_dp
        Quad%xi(7) = -0.2998304689007632_dp
        Quad%xi(8) = -0.1013262735219491_dp
        Quad%xi(9) = -Quad%xi(8)
        Quad%xi(10) = -Quad%xi(7)
        Quad%xi(11) = -Quad%xi(6)
        Quad%xi(12) = -Quad%xi(5)
        Quad%xi(13) = -Quad%xi(4)
        Quad%xi(14) = -Quad%xi(3)
        Quad%xi(15) = -Quad%xi(2)

        Quad%weights(2) = 0.05085036100592187_dp
        Quad%weights(3) = 0.08939369732593062_dp
        Quad%weights(4) = 0.1242553821325135_dp
        Quad%weights(5) = 0.15402698080716518_dp
        Quad%weights(6) = 0.17749191339170411_dp
        Quad%weights(7) = 0.19369002382520362_dp
        Quad%weights(8) = 0.20195830817822985_dp
        Quad%weights(9) = 0.20195830817822985_dp
        Quad%weights(10) = Quad%weights(7)
        Quad%weights(11) = Quad%weights(6)
        Quad%weights(12) = Quad%weights(5)
        Quad%weights(13) = Quad%weights(4)
        Quad%weights(14) = Quad%weights(3)
        Quad%weights(15) = Quad%weights(2)

    case(16)
        Quad%xi(2) = -0.9731321766314183_dp
        Quad%xi(3) = -0.910879995915574_dp
        Quad%xi(4) = -0.8156962512217703_dp
        Quad%xi(5) = -0.6910289806276847_dp
        Quad%xi(6) = -0.541385399330102_dp
        Quad%xi(7) = -0.3721744335654772_dp
        Quad%xi(8) = -0.189511973518317_dp
        Quad%xi(9) = 0.0_dp
        Quad%xi(10) = -Quad%xi(8)
        Quad%xi(11) = -Quad%xi(7)
        Quad%xi(12) = -Quad%xi(6)
        Quad%xi(13) = -Quad%xi(5)
        Quad%xi(14) = -Quad%xi(4)
        Quad%xi(15) = -Quad%xi(3)
        Quad%xi(16) = -Quad%xi(2)

        Quad%weights(2) = 0.04492194054325292_dp
        Quad%weights(3) = 0.07919827050368623_dp
        Quad%weights(4) = 0.11059290900702798_dp
        Quad%weights(5) = 0.13798774620192722_dp
        Quad%weights(6) = 0.1603946619976215_dp
        Quad%weights(7) = 0.1770042535156577_dp
        Quad%weights(8) = 0.18721633967761928_dp
        Quad%weights(9) = 0.19066187475346943_dp
        Quad%weights(10) = Quad%weights(8)
        Quad%weights(11) = Quad%weights(7)
        Quad%weights(12) = Quad%weights(6)
        Quad%weights(13) = Quad%weights(5)
        Quad%weights(14) = Quad%weights(4)
        Quad%weights(15) = Quad%weights(3)
        Quad%weights(16) = Quad%weights(2)

    case(17)
        Quad%xi(2) = -0.976105557412198_dp
        Quad%xi(3) = -0.920649185347533_dp
        Quad%xi(4) = -0.835593535218090_dp
        Quad%xi(5) = -0.723679329283243_dp
        Quad%xi(6) = -0.588504834318661_dp
        Quad%xi(7) = -0.434415036912123_dp
        Quad%xi(8) = -0.2663626528782805_dp
        Quad%xi(9) = -0.089749093484652_dp
        Quad%xi(10) = -Quad%xi(9)
        Quad%xi(11) = -Quad%xi(8)
        Quad%xi(12) = -Quad%xi(7)
        Quad%xi(13) = -Quad%xi(6)
        Quad%xi(14) = -Quad%xi(5)
        Quad%xi(15) = -Quad%xi(4)
        Quad%xi(16) = -Quad%xi(3)
        Quad%xi(17) = -Quad%xi(2)

        Quad%weights(2) = 0.039970628810914184_dp
        Quad%weights(3) = 0.07063716688563393_dp
        Quad%weights(4) = 0.0990162717175025_dp
        Quad%weights(5) = 0.12421053313296582_dp
        Quad%weights(6) = 0.1454119615738022_dp
        Quad%weights(7) = 0.16193951723760272_dp
        Quad%weights(8) = 0.17326210948945625_dp
        Quad%weights(9) = 0.17901586343970305_dp
        Quad%weights(10) = 0.17901586343970305_dp
        Quad%weights(11) = Quad%weights(8)
        Quad%weights(12) = Quad%weights(7)
        Quad%weights(13) = Quad%weights(6)
        Quad%weights(14) = Quad%weights(5)
        Quad%weights(15) = Quad%weights(4)
        Quad%weights(16) = Quad%weights(3)
        Quad%weights(17) = Quad%weights(2)

    case(18)
        Quad%xi(2) = -0.978611766222080_dp
        Quad%xi(3) = -0.928901528152586_dp
        Quad%xi(4) = -0.852460577796646_dp
        Quad%xi(5) = -0.751494202552613_dp
        Quad%xi(6) = -0.628908137265221_dp
        Quad%xi(7) = -0.488229285680714_dp
        Quad%xi(8) = -0.333504847824499_dp
        Quad%xi(9) = -0.169186023409282_dp
        Quad%xi(10) = 0.0_dp
        Quad%xi(11) = -Quad%xi(9)
        Quad%xi(12) = -Quad%xi(8)
        Quad%xi(13) = -Quad%xi(7)
        Quad%xi(14) = -Quad%xi(6)
        Quad%xi(15) = -Quad%xi(5)
        Quad%xi(16) = -Quad%xi(4)
        Quad%xi(17) = -Quad%xi(3)
        Quad%xi(18) = -Quad%xi(2)

        Quad%weights(2) = 0.035793365186175874_dp
        Quad%weights(3) = 0.0633818917626272_dp
        Quad%weights(4) = 0.08913175709920798_dp
        Quad%weights(5) = 0.11231534147730572_dp
        Quad%weights(6) = 0.1322672804487499_dp
        Quad%weights(7) = 0.14841394259593893_dp
        Quad%weights(8) = 0.16029092404406128_dp
        Quad%weights(9) = 0.16755658452714284_dp
        Quad%weights(10) = 0.17000191928482725_dp
        Quad%weights(11) = Quad%weights(9)
        Quad%weights(12) = Quad%weights(8)
        Quad%weights(13) = Quad%weights(7)
        Quad%weights(14) = Quad%weights(6)
        Quad%weights(15) = Quad%weights(5)
        Quad%weights(16) = Quad%weights(4)
        Quad%weights(17) = Quad%weights(3)
        Quad%weights(18) = Quad%weights(2)

    case(19)
        Quad%xi(2) = -0.980743704893914_dp
        Quad%xi(3) = -0.935934498812665_dp
        Quad%xi(4) = -0.866877978089950_dp
        Quad%xi(5) = -0.775368260952056_dp
        Quad%xi(6) = -0.663776402290311_dp
        Quad%xi(7) = -0.534992864031886_dp
        Quad%xi(8) = -0.392353183713909_dp
        Quad%xi(9) = -0.239551705922986_dp
        Quad%xi(10) = -0.080545937238822_dp
        Quad%xi(11) = -Quad%xi(10)
        Quad%xi(12) = -Quad%xi(9)
        Quad%xi(13) = -Quad%xi(8)
        Quad%xi(14) = -Quad%xi(7)
        Quad%xi(15) = -Quad%xi(6)
        Quad%xi(16) = -Quad%xi(5)
        Quad%xi(17) = -Quad%xi(4)
        Quad%xi(18) = -Quad%xi(3)
        Quad%xi(19) = -Quad%xi(2)

        Quad%weights(2) = 0.03223712318848816_dp
        Quad%weights(3) = 0.05718180212756649_dp
        Quad%weights(4) = 0.08063176399612001_dp
        Quad%weights(5) = 0.10199149969945108_dp
        Quad%weights(6) = 0.12070922762867593_dp
        Quad%weights(7) = 0.1363004823587244_dp
        Quad%weights(8) = 0.1483615540709169_dp
        Quad%weights(9) = 0.15658010264747546_dp
        Quad%weights(10) = 0.16074328638784577_dp
        Quad%weights(11) = 0.16074328638784577_dp
        Quad%weights(12) = Quad%weights(9)
        Quad%weights(13) = Quad%weights(8)
        Quad%weights(14) = Quad%weights(7)
        Quad%weights(15) = Quad%weights(6)
        Quad%weights(16) = Quad%weights(5)
        Quad%weights(17) = Quad%weights(4)
        Quad%weights(18) = Quad%weights(3)
        Quad%weights(19) = Quad%weights(2)

    case(20)
        Quad%xi(2) = -0.982572296604548_dp
        Quad%xi(3) = -0.941976296959746_dp
        Quad%xi(4) = -0.879294755323591_dp
        Quad%xi(5) = -0.796001926077712_dp
        Quad%xi(6) = -0.694051026062223_dp
        Quad%xi(7) = -0.575831960261831_dp
        Quad%xi(8) = -0.444115783279002_dp
        Quad%xi(9) = -0.301989856508765_dp
        Quad%xi(10) = -0.152785515802186_dp
        Quad%xi(11) = 0.0_dp
        Quad%xi(12) = -Quad%xi(10)
        Quad%xi(13) = -Quad%xi(9)
        Quad%xi(14) = -Quad%xi(8)
        Quad%xi(15) = -Quad%xi(7)
        Quad%xi(16) = -Quad%xi(6)
        Quad%xi(17) = -Quad%xi(5)
        Quad%xi(18) = -Quad%xi(4)
        Quad%xi(19) = -Quad%xi(3)
        Quad%xi(20) = -Quad%xi(2)

        Quad%weights(2) = 0.029184840098506866_dp
        Quad%weights(3) = 0.05184316900084789_dp
        Quad%weights(4) = 0.07327391818507369_dp
        Quad%weights(5) = 0.09298546795788497_dp
        Quad%weights(6) = 0.1105170832191237_dp
        Quad%weights(7) = 0.12545812119086924_dp
        Quad%weights(8) = 0.13745846286004137_dp
        Quad%weights(9) = 0.14623686244797748_dp
        Quad%weights(10) = 0.1515875751116814_dp
        Quad%weights(11) = 0.15338519033217496_dp
        Quad%weights(12) = Quad%weights(10)
        Quad%weights(13) = Quad%weights(9)
        Quad%weights(14) = Quad%weights(8)
        Quad%weights(15) = Quad%weights(7)
        Quad%weights(16) = Quad%weights(6)
        Quad%weights(17) = Quad%weights(5)
        Quad%weights(18) = Quad%weights(4)
        Quad%weights(19) = Quad%weights(3)
        Quad%weights(20) = Quad%weights(2)
        
    end select

end subroutine Spectral1DQuadrature

subroutine Spectral2DQuadrature(Quad, QuadBound)
    type(t_Quadrature)  :: Quad
    type(t_Quadrature)  :: QuadBound
    integer                             :: ii, jj, k
    
    Quad%n_points = QuadBound%n_points ** 2

    allocate(Quad%xi(Quad%n_points))
    allocate(Quad%eta(Quad%n_points))
    allocate(Quad%weights(Quad%n_points))

    k = 1
    do ii = 1, QuadBound%n_points
        do jj = 1, QuadBound%n_points
            Quad%xi(k) = QuadBound%xi(ii)
            Quad%eta(k) = QuadBound%xi(jj)
            Quad%weights(k) = QuadBound%weights(ii) * QuadBound%weights(jj)
            k = k + 1
        end do
    end do
    
end subroutine Spectral2DQuadrature

end module m_quadrature
