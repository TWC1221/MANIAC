program main
  use, intrinsic :: iso_fortran_env, only: output_unit
  use m_iter
  use m_PCG_solver, only: PCG_PRECON_NONE, PCG_PRECON_JACOBI, PCG_PRECON_CHOLESKY, PCG_PRECON_ILU, PETsc
  implicit none

  integer :: PCG_mode, max_iter
  integer :: t_start, t0, clock_rate
  real(8) :: cum_sec

  PCG_mode = PCG_PRECON_NONE
  max_iter = 200000

  ! Initialize cumulative timing and capture clock rate once
  call system_clock(t_start, clock_rate)
  cum_sec = 0.0d0

  ! Banner now includes solver name
  call print_banner('2026.2', solver_name(PCG_mode))

  call stage_begin(t0)
  call iter_xyz_multi(PCG_mode, max_iter, fixed =.true.)
  call stage_end('fixed', t0, t_start, clock_rate, cum_sec)

  call stage_begin(t0)
  call iter_xyz_multi(PCG_mode, max_iter, subcritical =.true.)
  call stage_end('subcritical', t0, t_start, clock_rate, cum_sec)

  call stage_begin(t0)
  call iter_xyz_multi(PCG_mode, max_iter, eigenvalue =.true.)
  call stage_end('eigenvalue', t0, t_start, clock_rate, cum_sec)

  call stage_begin(t0)
  call iter_xyz_multi(PCG_mode, max_iter, adjoint =.true.)
  call stage_end('adjoint', t0, t_start, clock_rate, cum_sec)

  call print_footer()

contains

  ! Map the integer PCG_mode to a readable solver name
  function solver_name(mode) result(name)
    integer, intent(in) :: mode
    character(len=32)   :: name
    select case (mode)
    case (PCG_PRECON_NONE);   name = 'PCG_PRECON_NONE'
    case (PCG_PRECON_JACOBI); name = 'PCG_PRECON_JACOBI'
    case (PCG_PRECON_CHOLESKY);   name = 'PCG_PRECON_CHOLESKY'
    case (PCG_PRECON_ILU);    name = 'PCG_PRECON_ILU'
    case (PETsc);             name = 'PETSc Wrapped PCG'
    case default;             name = 'UNKNOWN_SOLVER'
    end select
  end function solver_name

  subroutine print_banner(version, solver)
    use, intrinsic :: iso_fortran_env, only: output_unit
    character(len=*), intent(in) :: version, solver
    character(len=8)  :: date, time
    character(len=7)  :: tag_ddmonyy  ! e.g., 07Jan26
    character(len=8)  :: t_hms        ! hh:mm:ss

    call date_and_time(date=date, time=time)
    tag_ddmonyy = ddmonyy(date)
    t_hms = time(1:2)//':'//time(3:4)//':'//time(5:6)

    ! Left: midas, version, date-tag, solver; Right: current time
    write(output_unit,'(1X,"midas ",A,2X, A,2X,"[solver: ",A,"]", T73, A)') &
         trim(version), tag_ddmonyy, trim(solver), t_hms
    write(output_unit,'(1X,79("*"))')
    write(output_unit,*)
    call flush(output_unit)
  end subroutine

  pure function ddmonyy(date) result(tag)
    character(len=8), intent(in) :: date       ! "YYYYMMDD"
    character(len=7) :: tag                    ! "DDMonYY"
    character(len=2) :: dd, yy
    character(len=3) :: mon
    integer :: m

    dd = date(7:8)
    yy = date(3:4)
    read(date(5:6), '(I2)') m
    select case (m)
    case (1);  mon='Jan'
    case (2);  mon='Feb'
    case (3);  mon='Mar'
    case (4);  mon='Apr'
    case (5);  mon='May'
    case (6);  mon='Jun'
    case (7);  mon='Jul'
    case (8);  mon='Aug'
    case (9);  mon='Sep'
    case (10); mon='Oct'
    case (11); mon='Nov'
    case (12); mon='Dec'
    end select
    tag = dd//mon//yy
  end function

  subroutine print_footer()
    write(output_unit,'(1X,77("*"))')
    call flush(output_unit)
  end subroutine

  subroutine stage_begin(c0)
    integer, intent(out) :: c0
    call system_clock(c0)
  end subroutine

  subroutine stage_end(name, c0, t_start, clock_rate, cum_sec)
    character(len=*), intent(in) :: name
    integer,          intent(in) :: c0, t_start, clock_rate
    real(8),          intent(inout) :: cum_sec
    integer :: c1
    real(8) :: dt
    call system_clock(c1)
    dt = real(c1 - c0, 8) / real(clock_rate, 8)
    cum_sec = cum_sec + dt
    write(output_unit,'(1X, A, T73, F6.1, "s")') trim(name)//"...", cum_sec
    call flush(output_unit)
  end subroutine

end program main
!rm -rf build
!mkdir build && cd build
!cmake ..
!make -j


! # Delayed_Chi vector
! 0.03500755 0.1806982 0.1725102 0.3867821 0.1585752 0.0664266

! # Beta
! 0.0002284 0.0011787 0.001125 0.002523 0.0010344 0.0004333