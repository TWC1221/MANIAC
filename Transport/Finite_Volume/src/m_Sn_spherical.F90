!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  1D spherical S_N discrete ordinates transport solver (steady state)  -!
!  - Conservative shell balance with spherical face areas and volumes   -!
!  - Gauss–Legendre quadrature (N-1 to avoid μ=0)                       -!
!  - Angular differencing coefficients (EE Lewis CMNT formulation)      -!
!  - Isotropic external source and isotropic scattering                 -!
!  - Source iteration with inward/outward diamond-difference sweeps     -!
!  - Vacuum outer boundary; center regularity via angular differencing  -!
!  - Outputs phi.dat (r, φ, ψ_n) and a gnuplot script                   -!
!  - Colors ψ_n curves by their μ with colorbar scaled to [-1, 1]       -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 14/01/26     T. Charlton    Initial 1D spherical S_N solver; file     -!
!                             output and μ-colored gnuplot visualization-!
!------------------------------------------------------------------------!

program sn_spherical
  use m_constants
  use m_quadrature
  implicit none
  type(t_Quadrature) :: Quad

  integer :: N, I, nn, ii
  integer :: iter, max_iter
  real(8) :: R

  real(8), allocatable :: sigs(:), sigt(:)

  real(8), allocatable :: rho_edge(:), area(:), volume(:)
  
  real(8), allocatable :: mu(:), w(:)
  real(8), allocatable :: alpha(:)

  real(8), allocatable :: psi(:,:), psi_mu(:,:), psi_rho(:,:), psi_dir(:)
  real(8), allocatable :: S(:,:), q(:,:), phi(:)

  integer :: udat, uscr, col
  character(len=256) :: cmd
  character(len=32)  :: colstr, mustr
  real(dp) :: r_center, muval

  call system('pkill gnuplot')

  R = 1.0
  N = 20
  I = 1000

  allocate(sigs(I), sigt(I))
  sigs = 0.9
  sigt = 1.0

  call GetLineQuad(Quad, N-1)
  allocate(mu(N), w(N))
  mu = Quad%Xi
  w  = Quad%w

  allocate(rho_edge(I+1))
  rho_edge = [((R*ii)/(I), ii = 0, I)]

  allocate(area(I+1))
  area = 4 * pi * rho_edge**2

  allocate(volume(I))
  volume = [(4.0 * pi / 3.0 * ((rho_edge(ii+1))**3-(rho_edge(ii))**3), ii = 1, I)]

  allocate(alpha(N+1))
  alpha(1) = 0.0
  do nn = 1,N
    alpha(nn+1) = alpha(nn) - mu(nn)*w(nn)
    if (alpha(nn+1) <= 10E-6) alpha(nn+1) = 0.0
  end do 

  allocate(phi(I), psi(N,I),psi_mu(N+1,I),psi_rho(N,I+1),psi_dir(I+1), q(N,I),S(N,I))
  phi     = 0.0
  psi     = 0.0
  psi_rho = 0.0
  psi_mu  = 0.0
  psi_dir = 0.0

  ! ==============================
  ! Iterative Loop
  ! ==============================
  max_iter = 20

  do iter = 1, max_iter

    S = 1.0
    q = 0.0
    do ii = 1, I
      do nn = 1, N
        q(nn,ii) = sigs(ii) * phi(ii) + S(nn, ii)
      end do
    end do

    psi_rho = 0.0
    psi_mu  = 0.0
    psi_dir = 0.0

    do ii = I, 1, -1
      psi_mu(1,ii)=(2.0*psi_dir(ii+1)+(rho_edge(ii+1)-rho_edge(ii))*q(1,ii))/(2.0 + sigt(ii)*(rho_edge(ii+1)-rho_edge(ii)))
      psi_dir(ii) = 2.0*psi_mu(1,ii) - psi_dir(ii+1)
    end do

    do ii = I, 1, -1
      do nn = 1, N/2
        psi(nn,ii) = (abs(mu(nn))*(area(ii+1)+area(ii))*psi_rho(nn,ii+1)+(1/w(nn))*(area(ii+1)-area(ii))*(alpha(nn+1)+alpha(nn))*psi_mu(nn,ii)+volume(ii)*q(nn,ii))/(2.0*abs(mu(nn))*area(ii)+(2.0/w(nn))*(area(ii+1)-area(ii))*alpha(nn+1)+volume(ii)*sigt(ii))
        psi_rho(nn,ii) = 2.0*psi(nn,ii) - psi_rho(nn,ii+1)
        psi_mu(nn+1,ii) = 2.0*psi(nn,ii) - psi_mu(nn,ii)
      end do
    end do

    do nn = 1,N/2
      psi_rho(N+1-nn,1) = psi_rho(nn,1)
    end do 
    
    do ii = 1, I
      do nn = N/2+1,N
        psi(nn,ii) = (abs(mu(nn))*(area(ii+1)+area(ii))*psi_rho(nn,ii)+(1/w(nn))*(area(ii+1)-area(ii))*(alpha(nn+1)+alpha(nn))*psi_mu(nn,ii)+volume(ii)*q(nn,ii))/(2.0*abs(mu(nn))*area(ii+1)+(2.0/w(nn))*(area(ii+1)-area(ii))*alpha(nn+1)+volume(ii)*sigt(ii))
        psi_rho(nn,ii+1) = 2.0*psi(nn,ii) - psi_rho(nn,ii)
        psi_mu(nn+1,ii) = 2.0*psi(nn,ii) - psi_mu(nn,ii)
      end do
    end do

    phi = 0.0
    do ii = 1, I
      do nn = 1, N
        phi(ii) = phi(ii) + 0.5 * w(nn) * psi(nn,ii)
      end do
    end do

  end do

  print*,phi(1), phi(I)

  ! ==============================
  ! Write phi and psi to file
  ! ==============================
  udat = 10
  open(unit=udat, file='phi.dat', status='replace', action='write')
  write(udat,'(A)') '# r_center  phi  psi_1... psi_N'

  do ii = 1, I
      r_center = 0.5_dp*(rho_edge(ii) + rho_edge(ii+1))
      write(udat,'(F12.6,1X,ES16.8)', advance='no') r_center, phi(ii)
      do nn = 1, N
          write(udat,'(1X,ES16.8)', advance='no') psi(nn, ii)
      end do
      write(udat,*)
  end do
  close(udat)

  ! ==============================
  ! Build Gnuplot script
  ! ==============================
  uscr = 20
  open(unit=uscr, file='plot_sn.gp', status='replace', action='write')

  write(uscr,'(A)') "set grid"
  write(uscr,'(A)') "set xlabel 'r'"
  write(uscr,'(A)') "set ylabel 'Flux'"
  write(uscr,'(A)') "unset key"                
  write(uscr,'(A)') "set term qt"
  write(uscr,'(A)') "set colorbox vertical"
  write(uscr,'(A)') "set cblabel 'mu'"
  write(uscr,'(A)') "set cbrange [-1:1]"
  write(uscr,'(A)') "set cbtics -1,0.5,1"
  write(uscr,'(A)') "set palette defined (0 'blue', 0.5 'grey', 1 'red')"

  write(uscr,'(A)', advance='no') "plot "

  do nn = 1, N
      col = 2 + nn                   ! psi column in phi.dat
      muval = mu(nn)                 ! actual μ for this direction
      write(colstr,'(I0)') col
      write(mustr,'(F12.8)') muval

      if (nn < N) then
          write(uscr,'(A)') "'phi.dat' using 1:"//trim(colstr)// &
              " with lines lw 1 lc palette cb "//trim(adjustl(mustr))//" notitle, \"
      else
          write(uscr,'(A)') "'phi.dat' using 1:"//trim(colstr)// &
              " with lines lw 1 lc palette cb "//trim(adjustl(mustr))//" notitle, \"
      end if
  end do

  ! Add phi (column 2) on top as thick black, no title
  write(uscr,'(A)') "     'phi.dat' using 1:2 with lines lw 3 lc rgb 'black' notitle"

  close(uscr)

  cmd = "gnuplot -persist plot_sn.gp"
  call execute_command_line(cmd)

end program sn_spherical