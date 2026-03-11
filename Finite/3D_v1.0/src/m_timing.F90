module m_timing
    use m_constants
    use m_types
    use omp_lib
    implicit none
    
    integer, parameter :: MAX_TIMERS = 50
    integer, parameter :: NAME_LEN = 32    
    type(t_timer), save :: timers(MAX_TIMERS)
    integer, save       :: num_timers = 0
    
    contains
    
    subroutine timer_start(label)
        character(len=*), intent(in) :: label
        integer :: i
        do i = 1, num_timers
            if (trim(timers(i)%name) == trim(label)) then
                if (.not. timers(i)%active) then
                    timers(i)%start_t = omp_get_wtime()
                    timers(i)%active = .true.
                    timers(i)%calls = timers(i)%calls + 1
                end if
                return
            end if
        end do
        if (num_timers < MAX_TIMERS) then
            num_timers = num_timers + 1
            timers(num_timers)%name = label
            timers(num_timers)%start_t = omp_get_wtime()
            timers(num_timers)%active = .true.
            timers(num_timers)%calls = 1
        end if
    end subroutine
    
    subroutine timer_stop(label)
        character(len=*), intent(in) :: label
        integer :: i
        real(dp) :: t_now
        t_now = omp_get_wtime()
        do i = 1, num_timers
            if (trim(timers(i)%name) == trim(label)) then
                if (timers(i)%active) then
                    timers(i)%accumulated_t = timers(i)%accumulated_t + (t_now - timers(i)%start_t)
                    timers(i)%active = .false.
                end if
                return
            end if
        end do
    end subroutine
    
    subroutine print_timing_report(verbose)
        logical, intent(in), optional :: verbose
        integer :: i
        real(dp) :: total_time, avg_time, pct
        character(len=NAME_LEN+6) :: display_name
        logical :: do_print
        
        do_print = .true.
        if (present(verbose)) do_print = verbose
        if (.not. do_print) return
        
        total_time = 0.0_dp
        do i = 1, num_timers
            if (trim(timers(i)%name) == 'Total Execution') then
                total_time = timers(i)%accumulated_t
                exit
            end if
        end do
        if (total_time < 1.0e-9_dp) total_time = 1.0_dp
        
        write(*,*)
        write(*,'(A)') " [ TIMING ] :: Performance Analysis (Wall Clock)"
        write(*,'(A)') " |==========================================================================================|"
        write(*,'(A2, T5, A, T34, A8, T44, A12, T61, A12, T80, A10, T93, A1)') &
            " |", "Section", "Calls", "Time (s)", "Avg (s)", "Total (%)", "|"
        write(*,'(A)') " |==========================================================================================|"
        
        do i = 1, num_timers
            avg_time = 0.0_dp
            if (timers(i)%calls > 0) avg_time = timers(i)%accumulated_t / real(timers(i)%calls, dp)
            pct = (timers(i)%accumulated_t / total_time) * 100.0_dp
            
            select case (trim(timers(i)%name))
            case ("PETSc Init", "Mesh Parsing", "FE Generation", "Material Init", "Boundary Init", &
                  "Source Evaluation", "Linear Solve", "Eigen Update")
                
                display_name = "  >> " // trim(timers(i)%name)
                
                write(*,'(A2, T5, A30, T33, I8, T43, F12.4, T60, F12.6, T78, A7, F6.2, A1, T93, A1)') &
                    " | ", display_name, timers(i)%calls, timers(i)%accumulated_t, avg_time, &
                    "   >>", pct, "%", "|"
            
            case default
                display_name = trim(timers(i)%name)
                
                write(*,'(A2, T5, A30, T33, I8, T43, F12.4, T60, F12.6, T79, F9.2, A1, T93, A1)') &
                    " | ", display_name, timers(i)%calls, timers(i)%accumulated_t, avg_time, &
                    pct, "%", "|"
            end select
            
            if (trim(timers(i)%name) == 'Total Execution') then
                write(*,'(A)') " |------------------------------------------------------------------------------------------|"
            end if
        end do
        
        write(*,'(A)') " |==========================================================================================|"
        write(*,*)
    end subroutine
end module m_timing