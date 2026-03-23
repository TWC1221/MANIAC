!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  1D slab S_N discrete ordinates transport solver                      -!
!                                                                       -!
! Features:                                                             -!
!  - Steady-state fixed-source transport solution                       -!
!  - Isotropic scattering and isotropic external source                 -!
!  - Vacuum boundary conditions at RHS slab boundary                    -!
!  - Reflective boundary condition at LHS slab boundary                 -!
!  - Gauss–Legendre line quadrature (even-order, no μ = 0)              -!
!  - Uniform spatial mesh                                               -!
!                                                                       -!
! Spatial Discretizations:                                              -!
!  - Step differencing                                                  -!
!  - Diamond difference                                                 -!
!  - Linear discontinuous                                               -!
!                                                                       -!
! Numerical Methods:                                                    -!
!  - Source iteration                                                   -!
!  - Supports fixed-source and eigenvalue (k_eff) formulations          -!
!                                                                       -!
! Output:                                                               -!
!  - Scalar flux profiles written to text files                         -!
!    * flux_step.txt    : Step differencing solution                    -!
!    * flux_DD.txt      : Diamond-difference solution                   -!
!    * flux_LD.txt      : Linear-discontinuous solution                 -!
!  - Automatic gnuplot comparison of all methods                        -!
!                                                                       -!
!------------------------------------------------------------------------!
! Record of revisions:                                                  -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 11/01/26     T. Charlton    Initial 1D slab S_N solver with           -!
!                             diamond-difference; automated gnuplot     -!
! 12/02/26     T. Charlton    Implemented step method                   -!
! 20/02/26     T. Charlton    Major code refactor folowing S_N XYZ      -!
! 30/02/26     T. Charlton    Integrated LD method                      -!
! 02/20/26     T. Charlton    Implemented Reed Cell Benchmark           -!
!------------------------------------------------------------------------!

program sn_slab
  use m_constants
  use m_quadrature
  use m_outpVTK
  
  logical :: diamond
  integer :: I, N
  real(8) :: L
  real(8), allocatable :: sig_t(:), sig_s(:), Q(:)

  N = 16
  L = 8.0        
  I = 80

  allocate(sig_t(I), sig_s(I), Q(I))
  Sig_s(1:5*I/8) = 0.0
  Sig_s(5*I/8+1:I) = 0.9
  
  Sig_t(1:I/4) = 50.0
  Sig_t(I/4+1:3*I/8) = 5.0
  Sig_t(3*I/8+1:5*I/8) = 0.0
  Sig_t(5*I/8+1:I) = 1.0

  Q(1:I/4) = 50.0
  Q(I/4+1:3*I/8) = 0.0
  Q(3*I/8+1:5*I/8) = 0.0
  Q(5*I/8+1:3*I/4) = 1.0
  Q(3*I/4+1:I) = 0.0

  call sn_slab_DD_step(I, N, L, sig_s, sig_t, Q, diamond = .true.)
  call sn_slab_DD_step(I, N, L, sig_s, sig_t, Q, diamond = .false.)
  call sn_slab_LD(I, N, L, sig_s, sig_t, Q)

  call execute_command_line( &
  "gnuplot -persist -e ""set grid; " // &
  "set xlabel 'length (cm)'; set ylabel 'flux distribution (cm^{-2}s^{-1})'; " // &
  "plot " // &
  "'flux_step.txt' using 1:2 with lines lw 2 title 'step', " // &
  "'flux_DD.txt'   using 1:2 with lines lw 2 title 'diamond difference', " // &
  "'flux_LD.txt'   using 1:2 with lines lw 2 title 'linear discontinuous'""")
  
contains
subroutine sn_slab_DD_step(I, N, X, sigs, sigt, S_ext, diamond)
  implicit none
  type(t_quadrature) :: Quad

  integer, intent(in) :: N, I
  real(8), intent(in) :: X
  logical, intent(in) :: diamond

  integer :: nn, ii

  real(8), allocatable :: x_edge(:), dx(:)

  real(8), allocatable :: mu(:), w(:)

  real(8), allocatable, intent(in) :: sigs(:), sigt(:), S_ext(:)
  
  real(8), allocatable :: nu_sigf(:)

  real(8), allocatable :: psi(:,:), psi_i(:,:)

  real(8), allocatable :: q(:,:), phi(:), phi_prime(:), denom, mode
  
  integer :: istart, iend, istep, ifor, irev, f

  real(8) :: k_eff, k_eff_prime

  call system('pkill gnuplot')

  allocate(x_edge(I+1), dx(I))
  x_edge = [((X/I) * ii, ii = 0,I)]

  dx = [(x_edge(ii+1)-x_edge(ii),ii=1,I)]

  call GetLineQuad(Quad, N-1)
  mu = Quad%Xi
  w  = Quad%w

  allocate(nu_sigf(I))
  nu_sigf = 0.0

  allocate(psi(N,I),psi_i(N,I+1), q(N,I), phi(I), phi_prime(I))
  q = 0.0
  phi = 0.0
  phi_prime = 1.0
  k_eff = 0.0
  k_eff_prime = 1.0

  do
    do nn = 1,N
      do ii = 1,I
        q(nn,ii) = (sigs(ii) + nu_sigf(ii)/k_eff_prime) * phi_prime(ii) + S_ext(ii)
      end do
    end do

    do nn = 1,N/2

      istep  = merge(-1, 1, mu(nn) <= 0)
      istart = merge(I, 1, mu(nn)  <= 0)
      iend   = merge(1, I, mu(nn)  <= 0)

      ifor = merge(1, 0, mu(nn) >=0)
      irev = merge(1, 0, mu(nn) <=0)

      psi_i(nn, merge(1, I+1, mu(nn)  >= 0)) = 0.0
      
      mode = merge(2.0, 1.0, diamond)

      do ii = istart,iend,istep
        denom = sigt(ii) + (mode*abs(mu(nn))/dx(ii))
        psi(nn,ii) = ((mode*abs(mu(nn))/dx(ii))*psi_i(nn,ii+irev) + q(nn,ii))/denom

        if (diamond) then
          ! Diamond Difference
          psi_i(nn,ii+ifor) = 2.0 * psi(nn,ii) - psi_i(nn,ii+irev)
        else
          ! Step
          psi_i(nn,ii+ifor) = psi(nn,ii)
        end if

        phi(ii) = phi(ii) + 0.5*w(nn)*psi(nn,ii)
      end do

    end do
    
    do nn = 1,N/2
      psi_i(N+1-nn,1) = psi_i(nn,1)
    end do

    do nn = N/2+1,N

      istep  = merge(-1, 1, mu(nn) <= 0)
      istart = merge(I, 1, mu(nn)  <= 0)
      iend   = merge(1, I, mu(nn)  <= 0)

      ifor = merge(1, 0, mu(nn) >=0)
      irev = merge(1, 0, mu(nn) <=0)
      
      mode = merge(2.0, 1.0, diamond)

      do ii = istart,iend,istep
        denom = sigt(ii) + (mode*abs(mu(nn))/dx(ii))
        psi(nn,ii) = ((mode*abs(mu(nn))/dx(ii))*psi_i(nn,ii+irev) + q(nn,ii))/denom

        if (diamond) then
          ! Diamond Difference
          psi_i(nn,ii+ifor) = 2.0 * psi(nn,ii) - psi_i(nn,ii+irev)
        else
          ! Step
          psi_i(nn,ii+ifor) = psi(nn,ii)
        end if

        phi(ii) = phi(ii) + 0.5*w(nn)*psi(nn,ii)
      end do

    end do

    if  (maxval(S_ext) > 0.0 .and. maxval(nu_sigf) == 0) then
      if ( maxval(abs(phi - phi_prime)) < 1d-8 ) exit !FIXED SOURCE
      ! print*,maxval(abs(phi - phi_prime))
    else if (maxval(S_ext) == 0.0 .and. maxval(nu_sigf) > 0) then
      if ( maxval(abs(phi - phi_prime)) < 1d-8 .and. abs(k_eff - k_eff_prime) < 1d-10 ) exit !EIGENVALUE
      k_eff = k_eff_prime * sum(phi * nu_sigf)/sum(phi_prime * nu_sigf)
      k_eff_prime = k_eff
      ! print*,k_eff
    else
      print*,"Incompatible Inputs - S_external and nu_sigf mutually exclusive"
    end if 

    if (maxval(S_ext) > 0.0 .and. maxval(nu_sigf) == 0) phi_prime = phi !FIXED SOURCE
    if (maxval(S_ext) == 0.0 .and. maxval(nu_sigf) > 0) phi_prime = phi / sum(phi) !EIGENVALUE

    phi = 0.0

  end do

  ! print*,maxval(phi),minval(phi)

  if (diamond) open (action='write', file='flux_DD.txt', newunit=f, status='replace')
  if (.not. diamond) open (action='write', file='flux_step.txt', newunit=f, status='replace')

    do ii = 1,I
      write (f,'(ES16.8,1X,ES16.8)') ((ii-1)*dx(1)*I/(I-1)), phi(ii)
    end do 

  close(f)

end subroutine

subroutine sn_slab_LD(I, N, L, Sig_s, Sig_t, Q)
  implicit none

  type(t_quadrature) :: Quad

  integer, intent(in) :: I, N
  real(8), intent(in) :: L
  real(8), intent(in) :: Sig_t(:), Sig_s(:), Q(:)

  integer :: nn, ii, it, maxit, unit
  real(dp) :: h, tol, rel, eps
  real(dp), allocatable :: mu(:), w(:)
  real(dp), allocatable :: psi(:,:,:)
  real(dp), allocatable :: phiL(:), phiR(:), phiL_new(:), phiR_new(:)
  real(dp), allocatable :: CL(:), CR(:)
  real(dp), allocatable :: x(:), xc(:), phi_center(:)
  real(dp) :: incoming, A11, A12, A21, A22, b1, b2, det, psiL, psiR
  real(dp) :: nL_old, nR_old, nL_new, nR_new
  real(dp), allocatable :: inflow_L(:)
  integer :: nn_mirror

  allocate(inflow_L(N))
  inflow_L = 0.0_dp
          
  h     = L / real(I, dp)          
  tol   = 1.0e-6_dp       
  eps   = 1.0e-14_dp        
  maxit = 500 

  call GetLineQuad(Quad, N-1)
  mu = Quad%Xi
  w  = Quad%w

  allocate(psi(N, I, 2))     ! psi(nn, ii, 1=left; 2=right)
  allocate(phiL(I), phiR(I), phiL_new(I), phiR_new(I))
  allocate(CL(I), CR(I))
  allocate(x(0:I), xc(I), phi_center(I))

  psi  = 0.0_dp
  phiL = 0.0_dp
  phiR = 0.0_dp

  do ii = 0, I
     x(ii) = real(ii, dp) * h
  end do
  do ii = 1, I
     xc(ii) = 0.5_dp * (x(ii-1) + x(ii))
  end do

  do it = 1, maxit

    do ii = 1, I
      CL(ii) = (Sig_s(ii)/2.0_dp) * ( phiL(ii)*(h/3.0_dp) + phiR(ii)*(h/6.0_dp) ) + (Q(ii)/2.0_dp) * (h/2.0_dp)
      CR(ii) = (Sig_s(ii)/2.0_dp) * ( phiL(ii)*(h/6.0_dp) + phiR(ii)*(h/3.0_dp) ) + (Q(ii)/2.0_dp) * (h/2.0_dp)
    end do

    do nn = 1, N
      if (mu(nn) < 0.0_dp) then
        incoming = 0.0_dp
        do ii = I, 1, -1
          A11 = -mu(nn)/2.0_dp + Sig_t(ii)*(h/3.0_dp)
          A12 =  mu(nn)/2.0_dp + Sig_t(ii)*(h/6.0_dp)
          A21 = -mu(nn)/2.0_dp + Sig_t(ii)*(h/6.0_dp)
          A22 = -mu(nn)/2.0_dp + Sig_t(ii)*(h/3.0_dp)
          b1  = CL(ii)
          b2  = CR(ii) - mu(nn)*incoming
          det = A11*A22 - A12*A21
          psiL = ( b1*A22 - b2*A12 ) / det
          psiR = ( A11*b2 - A21*b1 ) / det
          psi(nn, ii, 1) = psiL
          psi(nn, ii, 2) = psiR
          incoming = psiL
        end do
      end if
    end do

    do nn = 1, N
      if (mu(nn) > 0.0_dp) then
        nn_mirror = N + 1 - nn
        if (nn_mirror >= 1.and. nn_mirror <= N.and. mu(nn_mirror) < 0.0_dp) then
          inflow_L(nn) = psi(nn_mirror, 1, 1)
        else
          inflow_L(nn) = 0.0_dp
        end if
      end if
    end do

    ! 3) Sweep μ > 0 (left -> right) using reflective inflow
    do nn = 1, N
      if (mu(nn) > 0.0_dp) then
        incoming = inflow_L(nn)
        do ii = 1, I
          A11 =  mu(nn)/2.0_dp + Sig_t(ii)*(h/3.0_dp)
          A12 =  mu(nn)/2.0_dp + Sig_t(ii)*(h/6.0_dp)
          A21 = -mu(nn)/2.0_dp + Sig_t(ii)*(h/6.0_dp)
          A22 =  mu(nn)/2.0_dp + Sig_t(ii)*(h/3.0_dp)
          b1  = CL(ii) + mu(nn)*incoming
          b2  = CR(ii)
          det = A11*A22 - A12*A21
          psiL = ( b1*A22 - b2*A12 ) / det
          psiR = ( A11*b2 - A21*b1 ) / det
          psi(nn, ii, 1) = psiL
          psi(nn, ii, 2) = psiR
          incoming = psiR
        end do
      end if
    end do

     phiL_new = 0.0_dp
     phiR_new = 0.0_dp
     do ii = 1, I
        do nn = 1, N
           phiL_new(ii) = phiL_new(ii) + w(nn) * psi(nn, ii, 1)
           phiR_new(ii) = phiR_new(ii) + w(nn) * psi(nn, ii, 2)
        end do
     end do

     nL_old = 0.0_dp
     nR_old = 0.0_dp
     nL_new = 0.0_dp
     nR_new = 0.0_dp
     do ii = 1, I
        nL_old = nL_old + phiL(ii)**2
        nR_old = nR_old + phiR(ii)**2
        nL_new = nL_new + phiL_new(ii)**2
        nR_new = nR_new + phiR_new(ii)**2
     end do
     nL_old = sqrt(nL_old)
     nR_old = sqrt(nR_old)
     nL_new = sqrt(nL_new)
     nR_new = sqrt(nR_new)
     rel = max( sqrt( sum( (phiL_new - phiL)**2 ) ) / ( nL_new + eps ), &
                sqrt( sum( (phiR_new - phiR)**2 ) ) / ( nR_new + eps ) )

     phiL = phiL_new
     phiR = phiR_new

     if (rel < tol) then
        print *, 'Converged in ', it, ' iterations. rel = ', rel
        exit
     end if
     if (it == maxit) then
        print *, 'Did not converge in ', maxit, ' iterations. rel = ', rel
     end if
  end do

  do ii = 1, I
     phi_center(ii) = 0.5_dp * (phiL(ii) + phiR(ii))
  end do

  open(newunit=unit, file='flux_LD.txt', status='replace', action='write', form='formatted')
  write(unit,'(A)') '# xc  phi_center  phi_L  phi_R'
  do ii = 1, I
     write(unit,'(4(1X,ES14.6))') (2.0*(ii-1)*xc(1)*I/(I-1)), phi_center(ii), phiL(ii), phiR(ii)
  end do
  close(unit)

end subroutine sn_slab_LD

end program sn_slab