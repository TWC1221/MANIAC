!------------------------------------------------------------------------!
! Program: sn_RZ                                                        -!
!                                                                       -!
! 2D axisymmetric (R–Z) S_N discrete ordinates neutron transport solver -!
!                                                                       -!
! - Steady-state transport equation                                     -!
! - Finite-volume discretization in cylindrical coordinates             -!
! - Diamond-difference spatial differencing                             -!
! - Even-order S_N angular quadrature (μ, η, ζ)                         -!
! - Isotropic scattering and optional fixed or fission source           -!
! - Source iteration and power iteration for k_eff                      -!
! - Vacuum boundary conditions on all external surfaces                 -!
!                                                                       -!
! Outputs:                                                              -!
! - Cell-centered scalar flux φ(r,z)                                    -!
! - Full angular flux ψ_{p,q}(r,z)                                      -!
! - VTK and ASCII output; automated gnuplot visualization               -!
!                                                                       -!
!------------------------------------------------------------------------!
! Revision history:                                                     -!
! 02/02/26  T. Charlton  Initial Working Version R–Z S_N solver         -!
!------------------------------------------------------------------------!

program sn_RZ
  use m_constants
  use m_quadrature
  use m_outpVTK
  implicit none
  type(t_sn_quadrature) :: sn_Quad

  real(8), allocatable :: R, Z, r_edge(:), z_edge(:), dr(:), dz(:)

  real(8), allocatable :: A(:,:), B(:), V(:,:)

  integer :: I, ii, J, jj, SN, P, pp, qq

  integer, allocatable :: Q(:)

  real(8), allocatable :: mu(:,:), zeta(:), w(:,:)

  real(8), allocatable :: alpha(:,:)

  real(8), allocatable :: psi(:,:,:,:), psi_q(:,:,:,:), psi_i(:,:,:,:), psi_j(:,:,:,:), psi_i_temp(:,:,:,:)

  real(8), allocatable :: Q_emission(:,:,:,:), S_external(:,:,:,:), sigs(:,:), sigt(:,:), nu_sigf(:,:), phi(:,:), phi_prime(:,:)

  real(8) :: denom, K_eff, K_eff_prime

  integer :: jstart, jend, jstep, jfor, jrev

  logical :: flag_adjoint = .false.

  call system('pkill gnuplot')

  SN = 16
  I = 300
  J = 300
  R = 1
  Z = 1

  allocate(r_edge(I+1), z_edge(J+1))
  r_edge = [((R*ii)/(I), ii = 0, I)]
  z_edge = [((Z*jj)/(J), jj = 0, J)]

  allocate(dr(I), dz(J))
  dr = [(r_edge(ii+1)-r_edge(ii),ii=1,I)]
  dz = [(z_edge(jj+1)-z_edge(jj),jj=1,J)]

  allocate(A(I+1,J), B(I), V(I,J))
  A = 2.0_dp*pi * matmul( reshape(r_edge(1:I+1), [I+1,1]), reshape(dz(1:J), [1,   J]) )
  B = pi * ( r_edge(2:I+1)**2 - r_edge(1:I)**2 )
  V = matmul( reshape(B, [I,1]), reshape(dz(1:J), [1,J]) )

  call GetRZAngleQuadrature(sn_Quad, SN, flag_adjoint)
  mu = sn_Quad%mu_pq(:,:)
  zeta = sn_Quad%zeta_p(:)
  w = sn_Quad%w_pq(:,:)
  P = sn_Quad%P
  Q = sn_Quad%qq_len(:)

  allocate(alpha(P,maxval(Q)+1))
  alpha(:,1) = 0.0
  do pp = 1,P
    do qq = 1,Q(pp)
      alpha(pp,qq+1) = alpha(pp,qq) - 2.0 * mu(pp,qq)*w(pp,qq)
      if (alpha(pp,qq+1) <= 10E-6) alpha(pp,qq+1) = 0.0
    end do
  end do 

  allocate(sigs(I,J), sigt(I,J), nu_sigf(I,J))
  sigs = 0.9
  sigt = 1.0
  nu_sigf = 0.4
  
  allocate(psi(P,maxval(Q),I,J), psi_q(P,maxval(Q)+1,I,J), psi_i(P,maxval(Q),I+1,J), psi_j(P,maxval(Q),I,J+1), psi_i_temp(P,maxval(Q),I+1,J), Q_emission(P,maxval(Q),I,J), S_external(P,maxval(Q), I,J), phi(I,J), phi_prime(I,J))
  Q_emission = 0.0
  S_external = 0.0

  psi =        0.0
  psi_q =      0.0
  psi_i =      0.0
  psi_j =      0.0

  phi =        0.0
  phi_prime =  1.0

  k_eff =      0.0
  k_eff_prime= 1.0

  do
    do pp = 1,P 
      do qq = 1,Q(pp)
        do ii = 1,I
          do jj = 1,J
            Q_emission(pp,qq,ii,jj) = (sigs(ii,jj) + nu_sigf(ii,jj)/k_eff_prime) * phi_prime(ii,jj) + S_external(pp,qq,ii,jj)
          end do
        end do
      end do
    end do

    do pp = 1,P ! psi_q(q=1/2) 

      jstep  = merge(-1, 1, zeta(pp) <= 0)
      jstart = merge(J, 1, zeta(pp)  <= 0)
      jend   = merge(1, J, zeta(pp)  <= 0)

      jfor = merge(1, 0, zeta(pp) >=0)
      jrev = merge(1, 0, zeta(pp) <=0)

      do ii = I, 1, -1
        do jj = jstart, jend, jstep
          denom = 2.0*((1.0-(zeta(pp)**2))**(0.5))/dr(ii) + 2.0*abs(zeta(pp))/dz(jj) + sigt(ii,jj)
          psi_q(pp,1,ii,jj) = (Q_emission(pp,1,ii,jj) + 2.0*((1.0-(zeta(pp)**2))**(0.5))*psi_i(pp,1,ii+1,jj)/dr(ii) + 2.0*abs(zeta(pp))*psi_j(pp,1,ii,jj+jrev)/dz(jj))/ denom

          psi_i(pp,1,ii,jj) = 2.0 * psi_q(pp,1,ii,jj) - psi_i(pp,1,ii+1,jj)
          psi_j(pp,1,ii,jj+jfor) = 2.0 * psi_q(pp,1,ii,jj) - psi_j(pp,1,ii,jj+jrev)
        end do
      end do

    end do 

    do pp = 1,P
      do qq = 1,Q(pp)/2 ! mu < 0

        jstep  = merge(-1, 1, zeta(pp) <= 0)
        jstart = merge(J, 1, zeta(pp)  <= 0)
        jend   = merge(1, J, zeta(pp)  <= 0)

        jfor = merge(1, 0, zeta(pp) >=0)
        jrev = merge(1, 0, zeta(pp) <=0)

        psi_i(pp, qq, I+1, :) = 0.0
        psi_j(pp, qq, :, J+1) = 0.0
        psi_j(pp, qq, :, 1) = 0.0

        do ii = I, 1, -1
          do jj = jstart, jend, jstep
            denom = (abs(mu(pp,qq))*(A(ii+1,jj)+A(ii,jj)) + 2.0*abs(zeta(pp))*B(ii) + (A(ii+1,jj)-A(ii,jj))*(alpha(pp,qq+1)+alpha(pp,qq))/(2.0*w(pp,qq)) + sigt(ii,jj)*V(ii,jj))
            psi(pp,qq,ii,jj) = (abs(mu(pp,qq))*(A(ii+1,jj)+A(ii,jj))*psi_i(pp,qq,ii+1,jj) + 2.0*abs(zeta(pp))*B(ii)*psi_j(pp,qq,ii,jj+jrev) + (A(ii+1,jj)-A(ii,jj))*(alpha(pp,qq+1)+alpha(pp,qq))*psi_q(pp,qq,ii,jj)/(2.0*w(pp,qq)) + V(ii,jj)*Q_emission(pp,qq,ii,jj)) / denom

            psi_i(pp,qq,ii,jj) = 2.0 * psi(pp,qq,ii,jj) - psi_i(pp,qq,ii+1,jj)
            psi_j(pp,qq,ii,jj+jfor) = 2.0 * psi(pp,qq,ii,jj) - psi_j(pp,qq,ii,jj+jrev)
            psi_q(pp,qq+1,ii,jj) = 2.0 * psi(pp,qq,ii,jj) - psi_q(pp,qq,ii,jj)

            phi(ii,jj) = phi(ii,jj) + 0.25*w(pp,qq)*psi(pp,qq,ii,jj)
          end do
        end do
      end do
    end do

    psi_i_temp = psi_i(:,:,:,:)
    do pp = 1, P ! Centerflux averaging for mu > 0
      do jj = 1, J
        where (mu(pp,1:Q(pp)) >= 0.0_dp) psi_i(pp,1:Q(pp),1,jj) = merge(sum( 2.0 * w(pp,1:Q(pp)) * abs(mu(pp,1:Q(pp))) * psi_i_temp(pp,1:Q(pp),1,jj), &
                    mask = mu(pp,1:Q(pp)) < 0.0_dp ) / sum( 2.0 * w(pp,1:Q(pp)) * mu(pp,1:Q(pp)), mask = mu(pp,1:Q(pp)) > 0.0_dp ), 0.0_dp, sum( 2.0 * w(pp,1:Q(pp)) * mu(pp,1:Q(pp)), mask = mu(pp,1:Q(pp)) > 0.0_dp ) > 0.0_dp)
      end do
    end do

    do pp = 1,P
      do qq = Q(pp)/2+1, Q(pp) ! mu > 0

        jstep  = merge(-1, 1, zeta(pp) <= 0)
        jstart = merge(J, 1, zeta(pp)  <= 0)
        jend   = merge(1, J, zeta(pp)  <= 0)

        jfor = merge(1, 0, zeta(pp) >=0)
        jrev = merge(1, 0, zeta(pp) <=0)

        psi_j(pp, qq, :, J+1) = 0.0
        psi_j(pp, qq, :, 1) = 0.0

        do ii = 1, I
          do jj = jstart, jend, jstep
            denom = (abs(mu(pp,qq))*(A(ii+1,jj)+A(ii,jj)) + 2.0*abs(zeta(pp))*B(ii) + (A(ii+1,jj)-A(ii,jj))*(alpha(pp,qq+1)+alpha(pp,qq))/(2.0*w(pp,qq)) + sigt(ii,jj)*V(ii,jj))
            psi(pp,qq,ii,jj) = (abs(mu(pp,qq))*(A(ii+1,jj)+A(ii,jj))*psi_i(pp,qq,ii,jj) + 2.0*abs(zeta(pp))*B(ii)*psi_j(pp,qq,ii,jj+jrev) + (A(ii+1,jj)-A(ii,jj))*(alpha(pp,qq+1)+alpha(pp,qq))*psi_q(pp,qq,ii,jj)/(2.0*w(pp,qq)) + V(ii,jj)*Q_emission(pp,qq,ii,jj)) / denom

            psi_i(pp,qq,ii+1,jj) = 2.0 * psi(pp,qq,ii,jj) - psi_i(pp,qq,ii,jj)
            psi_j(pp,qq,ii,jj+jfor) = 2.0 * psi(pp,qq,ii,jj) - psi_j(pp,qq,ii,jj+jrev)
            psi_q(pp,qq+1,ii,jj) = 2.0 * psi(pp,qq,ii,jj) - psi_q(pp,qq,ii,jj)

            phi(ii,jj) = phi(ii,jj) + 0.25*w(pp,qq)*psi(pp,qq,ii,jj)
          end do
        end do
      end do
    end do

    if  (maxval(S_external) > 0.0 .and. maxval(nu_sigf) == 0) then
      if ( maxval(abs(phi - phi_prime)) < 1d-8 ) exit !FIXED SOURCE
      print*,maxval(abs(phi - phi_prime))
    else if (maxval(S_external) == 0.0 .and. maxval(nu_sigf) > 0) then
      if ( maxval(abs(phi - phi_prime)) < 1d-8 .and. abs(k_eff - k_eff_prime) < 1d-10 ) exit !EIGENVALUE
      k_eff = k_eff_prime * sum(phi * nu_sigf)/sum(phi_prime * nu_sigf)
      k_eff_prime = k_eff
      print*,k_eff
    else
      print*,"Incompatible Inputs - S_external and nu_sigf mutually exclusive"
    end if 

    if (maxval(S_external) > 0.0 .and. maxval(nu_sigf) == 0) phi_prime = phi !FIXED SOURCE
    if (maxval(S_external) == 0.0 .and. maxval(nu_sigf) > 0) phi_prime = phi / sum(phi) !EIGENVALUE

    phi = 0.0

  end do
  
! -----------------------------
! OUTPUTS
! -----------------------------

call outpVTK_rectilinear('phi_rz.vtk', r_edge, z_edge, phi, I, J)
call outpVTK_rectilinear_vector_pq('SN_psi_decomposition_RZ.vtk', r_edge, z_edge, psi, I, J, P, Q, mu, zeta)

open(unit=10, file='phi_rz.dat', status='replace', action='write')
  do jj = 1, J
    do ii = 1, I
      write(10,'(F14.7,1X,F14.7,1X,ES16.8)') 0.5_dp * ( r_edge(ii) + r_edge(ii+1) ), 0.5_dp * ( z_edge(jj) + z_edge(jj+1) ), phi(ii,jj)
    end do
    write(10,*)
  end do
close(10)

open(unit=10, file='angles.dat', status='replace', action='write')
  do ii = 1, SN*(SN+2)/2
    write(10,'(3f20.10)') sn_quad%Angles(ii,1), &
                          sn_quad%Angles(ii,2), &
                          sn_quad%Angles(ii,3)
  end do
close(10)

call execute_command_line( &
"gnuplot -persist -e ""set term wxt; " // &
"set grid; " // &
"set xlabel 'x'; set ylabel 'y'; set cblabel 'phi'; set colorbox vertical; " // &
"set xrange [" // trim(adjustl(to_str(r_edge(1)))) // ":" // trim(adjustl(to_str(r_edge(I+1)))) // "]; " // &
"set yrange [" // trim(adjustl(to_str(z_edge(1)))) // ":" // trim(adjustl(to_str(z_edge(J+1)))) // "]; " // &
"set size ratio -1; " // &
"set view map; set pm3d at b; " // &
"set palette defined (0 '#000004', 0.25 '#3b0f70', 0.5 '#8c2981', 0.75 '#de4968', 1 '#fcfdbf'); " // &
"splot 'phi_rz.dat' using 1:2:3 with pm3d notitle""")

call execute_command_line( &
"gnuplot -persist -e ""set term wxt; " // &
"set grid; " // &
"set xlabel 'μ'; " // &
"set ylabel 'η'; " // &
"set zlabel 'ζ'; " // &
"set view equal xyz; " // &
"splot 'angles.dat' using 1:2:3 with points pt 7 ps 1.2 notitle""")

contains
  function to_str(x) result(s)
    real(dp), intent(in) :: x
    character(len=32) :: s
    write(s,'(F14.7)') x
  end function to_str

end program sn_RZ





