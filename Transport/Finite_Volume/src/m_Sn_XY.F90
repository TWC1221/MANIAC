!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  2D Cartesian S_N discrete ordinates transport solver                 -!
!  - Steady-state fixed-source solution                                 -!
!  - Eigenvalue (k-effective) solution                                  -!
!  - Adjoint (reverse/importance) transport capability                  -!
!  - Conservative cell balance on a uniform x-y grid                    -!
!  - Gauss–Legendre quadrature in 2D (μ, η angles)                      -!
!  - Diamond-difference (or step) spatial discretization                -!
!  - Isotropic external source and isotropic scattering                 -!
!  - Source iteration with inward/outward sweeps over all angles        -!
!  - Vacuum boundary conditions on all outer boundaries                 -!
!  - Outputs phi_xy.dat (x, y, φ) and plot_phi.gp gnuplot script        -!
!  - Visualizes scalar flux as a pm3d heatmap                           -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 17/01/26     T. Charlton    Initial 2D Cartesian S_N solver; file     -!
!                             output and pm3d gnuplot visualization     -!
!------------------------------------------------------------------------!

program sn_XY
  use m_constants
  use m_quadrature
  use m_outpVTK
  implicit none
  type(t_sn_quadrature) :: sn_Quad

  integer :: N, I, J, nn, ii, jj, Nang

  real(8), allocatable :: X, Y, x_edge(:), y_edge(:), x_center(:), y_center(:), dx(:), dy(:)

  real(8), allocatable :: mu(:), eta(:), w(:)

  real(8), allocatable :: sigs(:,:), sigt(:,:), nu_sigf(:,:)

  real(8), allocatable :: psi(:,:,:), psi_i(:,:,:), psi_j(:,:,:)

  real(8), allocatable :: S_ext(:,:,:), q(:,:,:), phi(:,:), phi_prime(:,:), denom, mode
  
  integer :: istart, iend, istep, jstart, jend, jstep, ifor, jfor, irev, jrev, udat, uscr

  real(8) :: k_eff, k_eff_prime

  logical :: diamond = .false., flag_adjoint = .false.

  call system('pkill gnuplot')

  N = 16
  I = 100
  J = 100
  X = 1.0
  Y = 1.0

  allocate(x_edge(I+1), y_edge(J+1), dx(I), dy(J))
  x_edge = [((X/I) * ii, ii = 0,I)]
  y_edge = [((Y/J) * jj, jj = 0,J)]

  dx = [(x_edge(ii+1)-x_edge(ii),ii=1,I)]
  dy = [(y_edge(jj+1)-y_edge(jj),jj=1,J)]

  call Get2DAngleQuadrature(sn_Quad, N, flag_adjoint)
  Nang = sn_Quad%NoAngles 
  mu = sn_Quad%Angles(1:Nang,1)
  eta = sn_Quad%Angles(1:Nang,2)
  w = sn_Quad%w(1:Nang)

  allocate(sigs(I,J), sigt(I,J), nu_sigf(I,J))
  sigs = 0.9
  sigt = 1.0
  nu_sigf = 0.0

  allocate(psi(Nang,I,J),psi_i(Nang,I+1,J),psi_j(Nang,I,J+1), q(Nang,I,J), S_ext(Nang,I,J), phi(I,J), phi_prime(I,J))
  q = 0.0
  S_ext = 1.0 
  phi = 0.0
  phi_prime = 1.0
  k_eff = 0.0
  k_eff_prime = 1.0

  do
    do nn = 1,Nang
      do jj = 1,J
        do ii = 1,I
          q(nn,ii,jj) = (sigs(ii,jj) + nu_sigf(ii,jj)/k_eff_prime) * phi_prime(ii,jj) + S_ext(nn,ii,jj)
        end do
      end do
    end do

    do nn = 1,Nang

      istep  = merge(-1, 1, mu(nn) <= 0)
      istart = merge(I, 1, mu(nn)  <= 0)
      iend   = merge(1, I, mu(nn)  <= 0)

      jstep  = merge(-1, 1, eta(nn) <= 0)
      jstart = merge(J, 1, eta(nn)  <= 0)
      jend   = merge(1, J, eta(nn)  <= 0)

      ifor = merge(1, 0, mu(nn) >=0)
      jfor = merge(1, 0, eta(nn) >=0)
      irev = merge(1, 0, mu(nn) <=0)
      jrev = merge(1, 0, eta(nn) <=0)

      psi_i(nn, merge(1, I+1, mu(nn)  >= 0), :) = 0.0
      psi_j(nn, :, merge(1, J+1, eta(nn) >= 0)) = 0.0

      mode = merge(2.0, 1.0, diamond)

      do jj = jstart,jend,jstep
        do ii = istart,iend,istep
          denom = sigt(ii,jj) + (mode*abs(mu(nn))/dx(ii)) + (mode*abs(eta(nn))/dy(jj))
          psi(nn,ii,jj) = ((mode*abs(mu(nn))/dx(ii))*psi_i(nn,ii+irev,jj) + (mode*abs(eta(nn))/dy(jj))*psi_j(nn,ii,jj+jrev) + q(nn,ii,jj))/denom
          if (diamond) then
            ! Diamond Difference
            psi_i(nn,ii+ifor,jj) = 2.0 * psi(nn,ii,jj) - psi_i(nn,ii+irev,jj)
            psi_j(nn,ii,jj+jfor) = 2.0 * psi(nn,ii,jj) - psi_j(nn,ii,jj+jrev) 
          else
            ! Step
            psi_i(nn,ii+ifor,jj) = psi(nn,ii,jj)
            psi_j(nn,ii,jj+jfor) = psi(nn,ii,jj)
          end if

          phi(ii,jj) = phi(ii,jj) + 0.25*w(nn)*psi(nn,ii,jj)
        end do
      end do

    end do
    
    if ( maxval(abs(phi - phi_prime)) < 1d-8 ) exit !FIXED SOURCE
    !if ( maxval(abs(phi - phi_prime)) < 1d-8 .and. abs(k_eff - k_eff_prime) < 1d-10 ) exit !EIGENVALUE

    !k_eff = k_eff_prime * sum(phi * nu_sigf)/sum(phi_prime * nu_sigf)
    !k_eff_prime = k_eff

    phi_prime = phi !FIXED SOURCE
    !phi_prime = phi / sum(phi) !EIGENVALUE

    phi = 0.0

  end do

  print*,maxval(phi),minval(phi), k_eff

! -----------------------------
! Write phi(x,y) as grid for pm3d
! -----------------------------

allocate(x_center(I), y_center(J))
do ii = 1, I
  x_center(ii) = 0.5_dp * ( x_edge(ii) + x_edge(ii+1) )
end do
do jj = 1, J
  y_center(jj) = 0.5_dp * ( y_edge(jj) + y_edge(jj+1) )
end do

udat = 10
open(unit=udat, file='phi_xy.dat', status='replace', action='write')
! Columns: x  y  phi
do jj = 1, J
  do ii = 1, I
    write(udat,'(F14.7,1X,F14.7,1X,ES16.8)') x_center(ii), y_center(jj), phi(ii,jj)
  end do
  write(udat,*) 
end do
close(udat)

! -----------------------------
! Build gnuplot script
! -----------------------------
uscr = 20
open(unit=uscr, file='plot_phi.gp', status='replace', action='write')

write(uscr,'(A)') "set term qt"
write(uscr,'(A)') "set grid"
write(uscr,'(A)') "set xlabel 'x'"
write(uscr,'(A)') "set ylabel 'y'"
write(uscr,'(A)') "set cblabel 'phi'"
write(uscr,'(A)') "set colorbox vertical"
write(uscr,'(A)') "set xrange ["//trim(adjustl(to_str(x_edge(1))))//":"//trim(adjustl(to_str(x_edge(I+1))))//"]"
write(uscr,'(A)') "set yrange ["//trim(adjustl(to_str(y_edge(1))))//":"//trim(adjustl(to_str(y_edge(J+1))))//"]"
write(uscr,'(A)') "set size ratio -1"            ! equal scaling x:y

write(uscr,'(A)') "set view map"
write(uscr,'(A)') "set pm3d at b"
write(uscr,'(A)') "set palette defined (0 '#000004', 0.25 '#3b0f70', 0.5 '#8c2981', 0.75 '#de4968', 1 '#fcfdbf')"
write(uscr,'(A)') "splot 'phi_xy.dat' using 1:2:3 with pm3d notitle"

write(uscr,'(A)') "# set contour base"
write(uscr,'(A)') "# set cntrparam levels 10"
write(uscr,'(A)') "# replot 'phi_xy.dat' using 1:2:3 with lines lc rgb 'black' notitle"

write(uscr,'(A)') "# unset view; set view 60,30"
write(uscr,'(A)') "# unset pm3d"
write(uscr,'(A)') "# splot 'phi_xy.dat' using 1:2:3 with lines lw 1 notitle"

close(uscr)

call execute_command_line("gnuplot -persist plot_phi.gp")

print *, '----------------------------------------'
print *, ' SN Quadrature'
print *, ' NoAngles = ', sn_quad%NoAngles
print *, '----------------------------------------'
print *, '  ii        mu           eta          zeta          w'
print *, '----------------------------------------'

do ii = 1, sn_quad%NoAngles
    print '(i4,5f13.6)', ii, &
          sn_quad%Angles(ii,1), &
          sn_quad%Angles(ii,2), &
          sn_quad%Angles(ii,3), &
          sn_quad%w(ii)
end do

print *, '----------------------------------------'

contains
  function to_str(x) result(s)
    real(dp), intent(in) :: x
    character(len=32) :: s
    write(s,'(F14.7)') x
  end function to_str

end program sn_XY