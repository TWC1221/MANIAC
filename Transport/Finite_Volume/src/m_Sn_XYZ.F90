!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  3D Cartesian S_N discrete ordinates transport solver                 -!
!  - Steady-state fixed-source solution                                 -!
!  - Eigenvalue (k-effective) solution                                  -!
!  - Adjoint (reverse/importance) transport capability                  -!
!  - Conservative cell balance on a uniform x-y-z grid                  -!
!  - Gauss–Legendre quadrature in 3D (μ, η, ζ angles)                   -!
!  - Diamond-difference (or step) spatial discretization                -!
!  - Isotropic external source and isotropic scattering                 -!
!  - Source iteration with inward/outward sweeps over all angles        -!
!  - Vacuum boundary conditions on all outer boundaries                 -!
!  - Outputs VTK files for scalar flux and angular flux vector fields   -!
!------------------------------------------------------------------------!
! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 27/01/26     T. Charlton    Initial 3D Cartesian S_N solver; VTK      -!
!                             output for φ and ψ                        -!
!------------------------------------------------------------------------!

program sn_XYZ
  use m_constants
  use m_quadrature
  use m_outpVTK
  implicit none
  type(t_sn_quadrature) :: sn_Quad

  integer :: N, I, J, K, nn, ii, jj, kk, Nang

  real(8), allocatable :: X, Y, Z, x_edge(:), y_edge(:), z_edge(:), dx(:), dy(:), dz(:)
  
  real(8), allocatable :: mu(:), eta(:), zeta(:), w(:)

  real(8), allocatable :: sigs(:,:,:), sigt(:,:,:), nu_sigf(:,:,:)

  real(8), allocatable :: psi(:,:,:,:), psi_i(:,:,:,:), psi_j(:,:,:,:), psi_k(:,:,:,:)

  real(8), allocatable :: S_ext(:,:,:,:), q(:,:,:,:), phi(:,:,:), phi_prime(:,:,:), denom, mode
  
  integer :: istart, iend, istep, jstart, jend, jstep, kstart, kend, kstep, ifor, jfor, kfor, irev, jrev, krev

  real(8) :: k_eff, k_eff_prime

  logical :: diamond = .true., flag_adjoint = .false.

  call system('pkill gnuplot')

  N = 16
  I = 10
  J = 10
  K = 10
  X = 1.0
  Y = 1.0
  Z = 1.0

  allocate(x_edge(I+1), y_edge(J+1), z_edge(K+1), dx(I), dy(J), dz(K))
  x_edge = [((X/I) * ii, ii = 0,I)]
  y_edge = [((Y/J) * jj, jj = 0,J)]
  z_edge = [((Z/K) * kk, kk = 0,K)]

  dx = [(x_edge(ii+1)-x_edge(ii),ii=1,I)]
  dy = [(y_edge(jj+1)-y_edge(jj),jj=1,J)]
  dz = [(z_edge(kk+1)-z_edge(kk),kk=1,K)]

  call Get3DAngleQuadrature(sn_Quad, N, flag_adjoint)
  Nang = sn_Quad%NoAngles 
  mu = sn_Quad%Angles(1:Nang,1)
  eta = sn_Quad%Angles(1:Nang,2)
  zeta = sn_Quad%Angles(1:Nang,3)
  w = sn_Quad%w(1:Nang)

  allocate(sigs(I,J,K), sigt(I,J,K), nu_sigf(I,J,K))
  sigs = 0.9
  sigt = 1.0
  nu_sigf = 0.0

  allocate(psi(Nang,I,J,K),psi_i(Nang,I+1,J,K),psi_j(Nang,I,J+1,K), psi_k(Nang,I,J,K+1), q(Nang,I,J,K), S_ext(Nang,I,J,K), phi(I,J,K), phi_prime(I,J,K))
  q = 0.0
  S_ext = 1.0 
  phi = 0.0
  phi_prime = 1.0
  k_eff = 0.0
  k_eff_prime = 1.0

  do
    do nn = 1,Nang
      do kk = 1,K
        do jj = 1,J
          do ii = 1,I
            q(nn,ii,jj,kk) = (sigs(ii,jj,kk) + nu_sigf(ii,jj,kk)/k_eff_prime) * phi_prime(ii,jj,kk) + S_ext(nn,ii,jj,kk)
          end do
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

      kstep  = merge(-1, 1, zeta(nn) <= 0)
      kstart = merge(K, 1, zeta(nn)  <= 0)
      kend   = merge(1, K, zeta(nn)  <= 0)

      ifor = merge(1, 0, mu(nn) >=0)
      jfor = merge(1, 0, eta(nn) >=0)
      kfor = merge(1, 0, zeta(nn) >=0)
      irev = merge(1, 0, mu(nn) <=0)
      jrev = merge(1, 0, eta(nn) <=0)
      krev = merge(1, 0, zeta(nn) <=0)

      psi_i(nn, merge(1, I+1, mu(nn)  >= 0), :, :) = 0.0
      psi_j(nn, :, merge(1, J+1, eta(nn) >= 0), :) = 0.0
      psi_k(nn, :, :, merge(1, K+1, zeta(nn) >= 0)) = 0.0

      mode = merge(2.0, 1.0, diamond)

      do kk = kstart,kend,kstep
        do jj = jstart,jend,jstep
          do ii = istart,iend,istep
            denom = sigt(ii,jj,kk) + (mode*abs(mu(nn))/dx(ii)) + (mode*abs(eta(nn))/dy(jj)) + (mode*abs(zeta(nn))/dz(kk))
            psi(nn,ii,jj,kk) = ((mode*abs(mu(nn))/dx(ii))*psi_i(nn,ii+irev,jj,kk) + (mode*abs(eta(nn))/dy(jj))*psi_j(nn,ii,jj+jrev,kk) + (mode*abs(zeta(nn))/dz(kk))*psi_k(nn,ii,jj,kk+krev) + q(nn,ii,jj,kk))/denom

            if (diamond) then
              ! Diamond Difference
              psi_i(nn,ii+ifor,jj,kk) = 2.0 * psi(nn,ii,jj,kk) - psi_i(nn,ii+irev,jj,kk)
              psi_j(nn,ii,jj+jfor,kk) = 2.0 * psi(nn,ii,jj,kk) - psi_j(nn,ii,jj+jrev,kk) 
              psi_k(nn,ii,jj,kk+kfor) = 2.0 * psi(nn,ii,jj,kk) - psi_k(nn,ii,jj,kk+krev) 
            else
              ! Step
              psi_i(nn,ii+ifor,jj,kk) = psi(nn,ii,jj,kk)
              psi_j(nn,ii,jj+jfor,kk) = psi(nn,ii,jj,kk)
              psi_k(nn,ii,jj,kk+kfor) = psi(nn,ii,jj,kk)
            end if

            phi(ii,jj,kk) = phi(ii,jj,kk) + 0.125*w(nn)*psi(nn,ii,jj,kk)
          end do
        end do           
      end do
    end do
    
    if ( maxval(abs(phi - phi_prime)) < 1d-8 ) exit !FIXED SOURCE
    !print*,maxval(abs(phi - phi_prime))
    ! if ( maxval(abs(phi - phi_prime)) < 1d-8 .and. abs(k_eff - k_eff_prime) < 1d-10 ) exit !EIGENVALUE
    ! k_eff = k_eff_prime * sum(phi * nu_sigf)/sum(phi_prime * nu_sigf)
    ! k_eff_prime = k_eff
    ! print*,k_eff

    phi_prime = phi !FIXED SOURCE
    !phi_prime = phi / sum(phi) !EIGENVALUE

    phi = 0.0

  end do

  call outpVTK_xyz_transport('SN_phi', phi, I, J, K, dx, dy, dz)
  call outpVTK_xyz_vector_field('SN_psi_field', psi, I, J, K, dx, dy, dz, Nang, mu, eta, zeta)
  call outpVTK_xyz_vector_decomposition('SN_psi_decomposition', psi, I, J, K, dx, dy, dz, Nang, mu, eta, zeta)

  !print*,maxval(phi),minval(phi), k_eff

  open(unit=10, file='angles.dat', status='replace')
    do ii = 1, Nang
      write(10,'(3f20.10)') sn_quad%Angles(ii,1), & ! mu eta zeta
                            sn_quad%Angles(ii,2), &
                            sn_quad%Angles(ii,3)
    end do
  close(10)

  call execute_command_line( &
  "gnuplot -persist -e ""set term wxt; " // &
  "set grid; " // &
  "set xlabel 'Angles(1)'; " // &
  "set ylabel 'Angles(2)'; " // &
  "set zlabel 'Angles(3)'; " // &
  "set view equal xyz; " // &
  "splot 'angles.dat' using 1:2:3 with points pt 7 ps 1.2 notitle""")

end program sn_XYZ