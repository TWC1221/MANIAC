module m_subcritical
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use mesh_types
  use CSR_types
  implicit none
contains

!======================================================================
!  subcriticality calculation
!======================================================================
subroutine k_subcritical(mesh, MATs, K_eff)

  type(MeshGrid),  intent(in) :: mesh
  type(MatMatrix), intent(in) :: MATs
  real(dp),        intent(in) :: K_eff
  real(dp)                 :: ks, kq, ksq, b_eff
  real(dp)                 :: ks_num, kq_num, S, Q, check

  integer :: G, ncell, gg
  real(dp), allocatable :: phi(:,:), phi_adj(:,:), phi_fix(:,:), src_ext(:,:)
  real(dp), allocatable :: nu_sigma_f(:), chi(:)

  !------------------------------------------------------------
  ! Dimensions
  !------------------------------------------------------------
  G     = MATs%G
  ncell = mesh%nx * mesh%ny * mesh%nz
  nu_sigma_f = MATs%nu_sigma_f
  chi         = MATs%chi

  allocate(phi(G,ncell), phi_adj(G,ncell), phi_fix(G,ncell), src_ext(G,ncell))

  !------------------------------------------------------------
  ! Read fields
  !------------------------------------------------------------
  call read_vtk_scalars('../critical_flux.vtk',         G, ncell, phi)
  call read_vtk_scalars('../G_adjoint.vtk',             G, ncell, phi_adj)
  call read_vtk_scalars('../fixed_flux.vtk',            G, ncell, phi_fix)
  call read_vtk_scalars('../src_external.vtk',          G, ncell, src_ext)

  !------------------------------------------------------------
  ! Initialize accumulators
  !------------------------------------------------------------
  ks_num    = 0.0d0
  kq_num    = 0.0d0
  S         = 0.0d0
  Q         = 0.0d0
  check     = 0.0d0

  !------------------------------------------------------------
  ! Loop over energy groups and accumulate ks/kq and b_eff numerator/denominator
  ! b_eff_num uses beta only (no per-group delayed spectrum available)
  !------------------------------------------------------------
  do gg = 1,G
    ks_num = ks_num + sum( phi_adj(gg,:) * nu_sigma_f(gg) * phi_fix(gg,:) * mesh%dV )
    S      = S      + sum(                 nu_sigma_f(gg) * phi_fix(gg,:) * mesh%dV )

    kq_num = kq_num + sum( phi_adj(gg,:) * src_ext(gg,:) * mesh%dV )
    Q      = Q      + sum(                 src_ext(gg,:) * mesh%dV )

    check = sum( nu_sigma_f(gg) * phi_fix(gg,:) * mesh%dV )
  end do

  ! Keep current functionality (even though this overwrites ks_num)
  ks_num = sum(phi_adj(1,:) * nu_sigma_f(2) * phi_fix(2,:) * mesh%dV)

  ks   = ks_num / (S)
  kq   = kq_num / (Q)
  ksq  = (ks_num + kq_num) / (S + Q)
  check = check - ks_num - kq_num

  ! Effective delayed fraction: adjoint-weighted ratio; equals beta here
  !b_eff = sum( phi_adj(2,:) * sum(MATs%Beta) * sum(MATs%Delayed_Chi) * nu_sigma_f(2) * phi_fix(2,:) * mesh%dV ) / sum( phi_adj(2,:) * nu_sigma_f(2) * phi_fix(2,:) * mesh%dV )

  if (check > 1.0d-7) print'(A,F12.8)', 'WARNING: subcritical multiplication factor self-consistency NOT confirmed, residual error = ', check
  if (sum(MATs%EVAL) > 0.1) then
    print '(A,F12.8,A,F15.8)', ' k_sub_ks   = ', ks, '    : error (pcm)  =', (ks - MATs%EVAL(2))/MATs%EVAL(2)*1.0d5
    print '(A,F12.8,A,F15.8)', ' k_sub_kq   = ', kq, '    : error (pcm)  =', (kq - MATs%EVAL(3))/MATs%EVAL(3)*1.0d5
    print '(A,F12.8,A,F15.8)', ' k_sub_ksq  = ', ksq, '    : error (pcm)  =', (ksq - MATs%EVAL(4))/MATs%EVAL(4)*1.0d5
    print '(A,F12.8,A,F15.8)', ' S          = ', S,  '    : error (pcm)  =', (S - MATs%EVAL(5))/MATs%EVAL(5)*1.0d5
  else
    print '(A,F12.8)',         ' k_sub_ks   = ', ks
    print '(A,F12.8)',         ' k_sub_kq   = ', kq
    print '(A,F12.8)',         ' k_sub_ksq  = ', ksq
    print '(A,F12.8)',         ' S          = ', S
  end if
  !print '(A,F12.8)',         'b_eff        = ', b_eff

end subroutine k_subcritical

!======================================================================
!  Read multigroup scalar fields from VTK
!======================================================================
subroutine read_vtk_scalars(filename, G, ncell, phi)
  implicit none

  character(len=*), intent(in) :: filename
  integer,          intent(in) :: G, ncell
  real(dp),         intent(out):: phi(G,ncell)

  character(len=256) :: line, tag
  integer :: i, gg, unit
  logical :: found(G)

  found = .false.
  unit  = 21

  open(unit, file=filename, status='old', action='read')

  do
    read(unit,'(A)', end=100) line
    do gg = 1, G
      write(tag,'(A,I0)') 'SCALARS phi_g', gg
      if (.not. found(gg) .and. index(line, trim(tag)) > 0) then
        read(unit,'(A)')
        do i = 1, ncell
          read(unit,*) phi(gg,i)
        end do
        found(gg) = .true.
      end if
    end do
    if (all(found)) exit
  end do

100 continue
  close(unit)

end subroutine read_vtk_scalars

end module m_subcritical
