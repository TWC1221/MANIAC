module m_driver
!------------------------------------------------------------------------!
!  Purpose:                                                             -!
!  Contains the subroutines required to drive solvers                   -!
!------------------------------------------------------------------------!
!  Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 09/10/2025    T. Charlton      Original code
! 15/10/2025    T. Charlton      Power Iterative  Method                -!
!------------------------------------------------------------------------! 
    use m_constants
    use m_MMS_solver
    use m_vacuum_BC_solver
    use m_iterative_power_solver
    use m_multigroup_1D_diffusion
    use m_Sparse_mm
    contains
    subroutine MMS_driver()
        implicit none
        integer :: N, ios, ii
        integer :: index = 0
        real(8), parameter :: alpha = 0.0d0  ! 0.0 = vacuum, 1.0 = reflective
        real(8), allocatable :: a(:), b(:), c(:)    ! tridiagonal matrix
        real(8), allocatable :: phi(:)             ! eigenvector
        real, allocatable :: x(:), y(:)
        real(8) :: L2Err, sum_x, sum_y, sum_xy, sum_xx, slope
        
        !-----------------------------
        ! Discretization Unit / Nodal Iteration Unit
        !-----------------------------
        ii=0
        open(unit=993, file="error.dat", status="replace", action="write")
        do N = 5, 5010, 1
            allocate(a(N), b(N), c(N), phi(N))
            call build_matrix_A_MMS(a, b, c, N, phi, L2Err)
            write(993,'(F10.5,1(1X,E15.8))') log10(real(N)), log10(L2Err)
            deallocate(a, b, c, phi)
            ii = ii + 1
        end do
        close(993)

        allocate(x(ii),y(ii))
        open(unit=993, file="error.dat", status='old', action='read', iostat=ios)
        do index = 1, ii 
            read(993, *, iostat=ios) x(index), y(index)
            sum_x = sum_x + x(index)
            sum_y = sum_y + y(index)
            sum_xy = sum_xy + x(index)*y(index)
            sum_xx = sum_xx + x(index)*x(index)
        end do
        close(993)

        slope = (ii*sum_xy - sum_x*sum_y) / (ii*sum_xx - sum_x*sum_x)

        print *, 'Log-Log Linear Fit: dy/dx = ', slope

        !-----------------------------
        ! Plot Error vs N
        !-----------------------------
        open(unit=994, file="plot_error.gp", status="replace", action="write")
        write(994,*) "set term wxt 0 title 'Error vs N'"
        write(994,*) "set title 'Error vs N'"
        write(994,*) "set xlabel 'log10(N)'"
        write(994,*) "set ylabel 'log10(Error)'"
        write(994,*) "plot 'error.dat' using 1:2 with lines title 'Error'"
        close(994)
        call system("gnuplot -persist plot_error.gp")

        !-----------------------------
        ! Plot Flux vs x
        !-----------------------------
        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'MMS derived Neutron Flux over 1-D domain (numerical & exact)'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical', \"
        write(995,*) "     'flux.dat' using 1:3 with lines title 'exact'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")

    end subroutine MMS_driver

    subroutine vacuum_BC_driver()
        implicit none
        integer :: N, ios, ii, index
        real(8), parameter :: alpha = 0.0d0  ! 0.0 = vacuum, 1.0 = reflective
        real(8), allocatable :: a(:), b(:), c(:)    ! tridiagonal matrix
        real(8), allocatable :: phi(:), phi_analytical(:)             ! eigenvector
        real, allocatable :: x(:), y(:)
        real(8) :: L2Err, sum_x, sum_y, sum_xy, sum_xx, slope
        
        !-----------------------------
        ! Discretization Unit / Nodal Iteration Unit
        !-----------------------------
        open(unit=993, file="error.dat", status="replace", action="write")
        ii=0
        do N = 5, 40, 1
            allocate(a(N), b(N), c(N), phi(N), phi_analytical(N))
            call build_matrix_A_vacuum(a, b, c, N, alpha, phi, phi_analytical, L2Err)
            write(993,'(F10.5,1(1X,E15.8))') log10(real(N)), log10(L2Err)
            deallocate(a, b, c, phi, phi_analytical)
            ii = ii + 1
        end do
        close(993)

        allocate(x(ii),y(ii))
        open(unit=993, file="error.dat", status='old', action='read', iostat=ios)
        do index = 1, ii 
            read(993, *, iostat=ios) x(index), y(index)
            sum_x = sum_x + x(index)
            sum_y = sum_y + y(index)
            sum_xy = sum_xy + x(index)*y(index)
            sum_xx = sum_xx + x(index)*x(index)
        end do
        close(993)

        slope = (ii*sum_xy - sum_x*sum_y) / (ii*sum_xx - sum_x*sum_x)

        print *, 'Log-Log Linear Fit: dy/dx = ', slope

        !-----------------------------
        ! Plot Error vs N
        !-----------------------------
        open(unit=994, file="plot_error.gp", status="replace", action="write")
        write(994,*) "set term wxt 0 title 'Error vs N'"
        write(994,*) "set title 'Error vs N'"
        write(994,*) "set xlabel 'log10(N)'"
        write(994,*) "set ylabel 'log10(Error)'"
        write(994,*) "plot 'error.dat' using 1:2 with lines title 'Error'"
        close(994)
        call system("gnuplot -persist plot_error.gp")

        !-----------------------------
        ! Plot Flux vs x
        !-----------------------------
        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical', \"
        write(995,*) "     'flux.dat' using 1:3 with lines title 'exact'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")

    end subroutine Vacuum_BC_driver

    subroutine heterogeneous_driver_m1()
        implicit none
        integer :: N = 1000, N_interface = 500
        real(8) :: alpha = 0
        real(8), allocatable :: a(:), b(:), c(:), phi(:)

        allocate(a(N), b(N), c(N), phi(N))

        call build_matrix_A_heterogeneous_m1(a, b, c, N, N_interface, alpha, phi)

        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x for different alpha'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")
    end subroutine heterogeneous_driver_m1

    subroutine heterogeneous_driver_m2()
        implicit none
        integer :: N = 100, N_interface = 5
        real(8) :: alpha = 0.0, x_interface = 0.4
        real(8), allocatable :: a(:), b(:), c(:), phi(:)

        allocate(a(N), b(N), c(N), phi(N))

        call build_matrix_A_heterogeneous_m2(a, b, c, N, N_interface, x_interface, alpha, phi)

        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x for different alpha'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical', 'flux.dat' using 1:3 with lines title 'analytical'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")
    end subroutine heterogeneous_driver_m2

    subroutine iterative_power_driver()
        implicit none
        integer :: N = 1000, jj = 1
        real(8) :: alpha = 0.5, x_interface = 1, nu = 3, Norm
        real(8), allocatable :: a(:), b(:), c(:), phi(:), dx(:), Dif(:), Sigma_a(:), Sigma_f(:), S0(:), phi_old(:)
        real(8), dimension(1000) :: lambda

        allocate(dx(N), Dif(N), Sigma_a(N), Sigma_f(N), S0(N), phi(N), a(N), b(N), c(N), phi_old(N))
        call system("pkill gnuplot")

        !call random_number(phi(1:N)) 
        phi(1:N) = 1 
        lambda(1) = 1; !LAMBDA 0

        dx(1:N) = x_interface / real(N, 8)

        Dif(1:N) = 1.0/(3.0) 
        Sigma_a(1:N) = 2
        Sigma_f(1:N) = 0.16

        do
            S0 = (nu * Sigma_f * phi) / lambda(jj) !S(0)
            S0(1:N) = S0(1:N)

            jj = jj + 1

            phi_old = phi !PHI 0

            call build_matrix_A_iterative_power(a, b, c, N, alpha, phi, dx, Dif, Sigma_a, S0)
            print*,phi

            lambda(jj) = lambda(jj-1) * sum(nu*Sigma_f*phi*dx)/sum(nu*Sigma_f*phi_old*dx)

            Norm = nu*sum(phi*Sigma_f*dx)/(lambda(jj))
            phi = phi/Norm

            print*,"Iteration =", jj, "Keff =", lambda(jj), "norm_check =",nu*sum(phi*Sigma_f*dx)/(lambda(jj))

            if (abs((lambda(jj) - lambda(jj-1)) / lambda(jj-1)) < 1.0d-10) exit
            if (maxval(abs(phi - phi_old)) < 1.0d-8) exit
        end do

        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x for different alpha'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "set yrange [0:*]"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")
    end subroutine iterative_power_driver

!--------------------------------------------------------------------------------------------!
!  Purpose:                                                                                 -!
!  Contains inner/outer iteration on multigroup diffusion eigenvalue problem                -!
!           no-upscatter fixed fission source solver                                        -! 
!--------------------------------------------------------------------------------------------!
!  Record of revisions:                                                                     -!
!   Date       Programmer     Description of change                                         -!
!   ====       ==========     =====================                                         -!
! 20/10/2025    T. W Charlton     Original code                                             -!
! 21/10/2025    T. W Charlton     Added Fixed Source & Upscatter Capabilities               -!
! 22/10/2025    T. W Charlton     Refined & Implemented Upscatter Fission Eigenvalue Solver -!
!--------------------------------------------------------------------------------------------! 

    subroutine multigroup_diffusion_iter(mode, PCG_mode)
        use CSR_types, only: CSRMatrix
        implicit none
        type(CSRMatrix) :: A_all(3)
        integer, intent(in) :: mode, PCG_mode ! 1 = Fixed Source No Upscatter, 2 = Fixed Source Upscatter Iteration, 3 = Multigroup Iterative Fission No Upscatter, 4 = Multigroup Iterative Fission with Upscatter Iteration
        integer :: N = 1000, G = 3, gg, jj, ii
        real(8) :: alpha = 0.0, R_domain = 10
        real(8), allocatable :: phi(:,:), phi_prime(:,:), phi_ptr(:), Dif(:), Sigma_t(:), Sigma_a(:), nu_Sigma_f(:), Sigma_s(:,:), Sigma_sr(:,:), S_f(:), S_f_iter(:), K_eff(:), x(:), dx(:), dx_V(:), chi(:)
        real(8), allocatable :: Sigma_s_upscatter(:,:), Sigma_s_L(:,:), Sigma_s_U(:,:), Sigma_s_L_sum(:,:), Sigma_s_U_sum(:,:) ! Experimental

        allocate(phi(G,N), phi_prime(G,N), phi_ptr(N), Dif(G), Sigma_t(G), Sigma_a(G), nu_Sigma_f(G), Sigma_s(G,G), Sigma_sr(G,G), S_f(N), S_f_iter(N), K_eff(1000), x(N), dx(N-1), dx_V(N), chi(G))
        allocate(Sigma_s_upscatter(1:G,1:G), Sigma_s_L(1:G,1:G), Sigma_s_L_sum(1:G,1:N), Sigma_s_U(1:G,1:G), Sigma_s_U_sum(1:G,1:N))
        call system("pkill gnuplot")

        dx(1:N-1) = R_domain / real(N-1, 8) 
        dx_V(1:N-1) = 0.5*dx(1:N-1) ; dx_V(2:N) = dx_V(2:N) + 0.5 * dx(1:N-1)

        chi = [0.9, 0.1, 0.0]
        Sigma_a = [0.015, 0.04, 0.12]
        nu_Sigma_f = [0.02, 0.1, 0.35] !0.35]

        S_f = 1 ; K_eff(1) = 1; 
        select case (mode)

        case (1) ! FIXED SOURCE NO-UPSCATTER
            print *, "FIXED SOURCE NO-UPSCATTER MODE"
            Sigma_s = reshape([0.20d0, 0.00d0, 0.00d0, 0.05d0, 0.25d0, 0.00d0, 0.0d0, 0.07d0, 0.3d0], shape=[3,3])
            do gg = 1,G
                Sigma_sr(gg,1:G) = Sigma_s(gg,1:G); Sigma_sr(gg,gg) = 0.0d0
                Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s(gg,1:G))
                Dif(gg) = 1.0d0/(3.0d0*Sigma_t(gg))
                call build_CSR_matrix_multigroup(N, alpha, G, gg, dx, Dif(gg), Sigma_s(1:G,1:G), Sigma_a(gg), A_all(gg))
            end do
            do gg = 1,3
                S_f_iter(1:N) = chi(gg) * S_f(1:N) + Sigma_sr(1,gg)*phi(1,1:N) + Sigma_sr(2,gg)*phi(2,1:N)
                phi_ptr = phi(gg,1:N)
                call PCG_algorithm(A_all(gg)%AA, A_all(gg)%JA, A_all(gg)%IA, phi_ptr, S_f, PCG_mode, N, dx)
                phi(gg,1:N) = phi_ptr
            end do

        case (2) ! FIXED SOURCE UP-SCATTER ITERATION
            print *, "FIXED SOURCE UP-SCATTER ITERATION MODE"
            Sigma_s_upscatter = reshape([0.20d0, 0.10d0, 0.17d0, 0.05d0, 0.25d0, 0.10d0, 0.10d0, 0.07d0, 0.3d0], shape=[3,3])
            Sigma_s_L(:,:) = 0 ; Sigma_s_L(2:3,1:2) = Sigma_s_upscatter(2:3,1:2) ; Sigma_s_L(2,2) = 0 ;
            Sigma_s_U(:,:) = 0 ; Sigma_s_U(1:2,2:3) = Sigma_s_upscatter(1:2,2:3) ; Sigma_s_U(2,2) = 0 ;

            !S_f = 0 !Surface_Source only problems

            do gg = 1,G
                Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s_upscatter(gg,1:G))
                Dif(gg) = 1.0d0/(3.0d0*Sigma_t(gg))
                call build_CSR_matrix_multigroup(N, alpha, G, gg, dx, Dif(gg), Sigma_s_upscatter(1:G,1:G), Sigma_a(gg), A_all(gg))
            end do

            do 
                do gg = 1,3
                    do jj = 1,3
                        Sigma_s_L_sum(jj,1:N) = Sigma_s_L(jj,gg) * phi_prime(jj,1:N)
                        Sigma_s_U_sum(jj,1:N) = Sigma_s_U(jj,gg) * phi(jj,1:N)
                    end do

                    S_f_iter(1:N) = chi(gg) * S_f(1:N) + sum(Sigma_s_L_sum,dim=1) + sum(Sigma_s_U_sum,dim=1)

                    phi_ptr = phi(gg,1:N)
                    call PCG_algorithm(A_all(gg)%AA, A_all(gg)%JA, A_all(gg)%IA, phi_ptr, S_f, PCG_mode, N, dx)!build_matrix_multigroup
                    phi(gg,1:N) = phi_ptr
                end do

                print'(10F10.5)',phi(1,N/2),phi(2,N/2),phi(3,N/2)

                if (maxval(abs(phi(3,1:N) - phi_prime(3,1:N))) < 1.0d-8 .and. maxval(abs(phi(2,1:N) - phi_prime(2,1:N))) < 1.0d-8) exit

                phi_prime(3,1:N) = phi(3,1:N)
                phi_prime(2,1:N) = phi(2,1:N) 
            end do

        case(3) ! MULTIGROUP ITERATIVE FISSION SOLVER NO-UPSCATTER
            print *, "MULTIGROUP ITERATIVE FISSION SOLVER NO-UPSCATTER MODE"
            Sigma_s = reshape([0.20d0, 0.00d0, 0.00d0, 0.05d0, 0.25d0, 0.00d0, 0.0d0, 0.07d0, 0.3d0], shape=[3,3])
            do gg = 1,G
                Sigma_sr(gg,1:G) = Sigma_s(gg,1:G); Sigma_sr(gg,gg) = 0.0d0
                Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s(gg,1:G))
                Dif(gg) = 1.0d0/(3.0d0*Sigma_t(gg))
                call build_CSR_matrix_multigroup(N, alpha, G, gg, dx, Dif(gg), Sigma_s_upscatter(1:G,1:G), Sigma_a(gg), A_all(gg))
            end do

            do jj = 2,10
                S_f_iter = S_f
                do gg = 1,3
                    S_f(1:N) = chi(gg)/K_eff(jj-1) * S_f_iter(1:N) + Sigma_sr(1,gg)*phi(1,1:N) + Sigma_sr(2,gg)*phi(2,1:N)
                    phi_ptr = phi(gg,1:N)
                    call PCG_algorithm(A_all(gg)%AA, A_all(gg)%JA, A_all(gg)%IA, phi_ptr, S_f, PCG_mode, N, dx)
                    phi(gg,1:N) = phi_ptr
                end do

                S_f(1:N) = matmul(transpose(phi),nu_Sigma_f)
                K_eff(jj) = K_eff(jj-1)*sum(S_f*dx_V)/sum(S_f_iter*dx_V)
    
                if (abs((K_eff(jj) - K_eff(jj-1)) / K_eff(jj-1)) < 1.0d-12) exit
                if (maxval(abs(S_f - S_f_iter)) < 1.0d-8) exit
            end do
            
            print*,"Iteration =", jj-1, "Keff =", K_eff(jj-1)

        case(4, 5) ! MULTIGROUP ITERATIVE FISSION SOLVER with Up-Scatter Iteration
            print *, "MULTIGROUP ITERATIVE FISSION SOLVER WITH UP-SCATTER ITERATION MODE"
            Sigma_s_upscatter = reshape([0.20d0, 0.10d0, 0.17d0, 0.05d0, 0.25d0, 0.10d0, 0.10d0, 0.07d0, 0.3d0], shape=[3,3]) ! TEST DATA (UNREALISTIC UPSCATTER)
            !Sigma_s_upscatter = reshape([0.30d0, 0.20d0, 0.085d0,0.02d0, 0.25d0, 0.09d0, 0.001d0, 0.009d0, 0.17d0], shape=[3,3]) ! REALISTIC DATA

            if (mode == 5) Sigma_s_upscatter = transpose(Sigma_s_upscatter)

            Sigma_s_L(:,:) = 0 ; Sigma_s_L(2:3,1:2) = Sigma_s_upscatter(2:3,1:2) ; Sigma_s_L(2,2) = 0 ;
            Sigma_s_U(:,:) = 0 ; Sigma_s_U(1:2,2:3) = Sigma_s_upscatter(1:2,2:3) ; Sigma_s_U(2,2) = 0 ;

            do gg = 1,G
                Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s_upscatter(gg,1:G))
                Dif(gg) = 1.0d0/(3.0d0*Sigma_t(gg))
                call build_CSR_matrix_multigroup(N, alpha, G, gg, dx, Dif(gg), Sigma_s_upscatter(1:G,1:G), Sigma_a(gg), A_all(gg))
            end do

            do ii = 2,1000
                S_f_iter(1:N) = S_f(1:N)
                do 
                    do gg = 1,3
                        do jj = 1,3
                            Sigma_s_L_sum(jj,1:N) = Sigma_s_L(jj,gg) * phi_prime(jj,1:N)
                            Sigma_s_U_sum(jj,1:N) = Sigma_s_U(jj,gg) * phi(jj,1:N)
                        end do

                        if (mode == 4) S_f(1:N) = chi(gg)/K_eff(ii-1) * S_f_iter(1:N) + sum(Sigma_s_L_sum,dim=1) + sum(Sigma_s_U_sum,dim=1)
                        if (mode == 5) S_f(1:N) = nu_sigma_f(gg)/K_eff(ii-1) * S_f_iter(1:N) + sum(Sigma_s_L_sum,dim=1) + sum(Sigma_s_U_sum,dim=1)

                        phi_ptr = phi(gg,1:N)
                        call PCG_algorithm(A_all(gg)%AA, A_all(gg)%JA, A_all(gg)%IA, phi_ptr, S_f, PCG_mode, N, dx)
                        phi(gg,1:N) = phi_ptr
                    end do

                    !print'(10F10.5)',phi(1,N/2),phi(2,N/2),phi(3,N/2)
                    if (maxval(abs(phi(3,1:N) - phi_prime(3,1:N))) < 1.0d-7 .and. maxval(abs(phi(2,1:N) - phi_prime(2,1:N))) < 1.0d-7) exit
                    phi_prime(3,1:N) = phi(3,1:N)
                    phi_prime(2,1:N) = phi(2,1:N) 
                    exit ! Accelerates convergence by performing many more outer iterations and negleting inner iterations
                end do

                if (mode == 4) S_f(1:N) = matmul(transpose(phi),nu_Sigma_f)
                if (mode == 5) S_f(1:N) = matmul(transpose(phi),chi)
                K_eff(ii) = K_eff(ii-1)*sum(S_f*dx_V)/sum(S_f_iter*dx_V)
                if (mode == 4) print*,"Iteration =", ii-1, "Keff =", K_eff(ii-1)
                if (abs((K_eff(ii) - K_eff(ii-1)) / K_eff(ii-1)) < 1.0d-6 .and. maxval(abs(S_f - S_f_iter)) < 1.0d-6) exit
            end do

            !if (mode == 4) print*,"Iteration =", ii-1, "Keff =", K_eff(ii-1)

        case default
            print *, "ERROR: INVALID MODE SELECTED"
        end select

        x(1) = 0.0d0
        do jj = 2, N
            x(jj) =  sum(dx(1:jj-1))
        end do

        open(unit=991, file="flux.dat", status="replace", action="write")
        do jj = 1, N
            write(991,'(F10.5,3(1X,E15.8))') x(jj), phi(1,jj), phi(2,jj), phi(3,jj)
        end do
        close(991)

        open(unit=991, file="plot_flux.gp", status="replace", action="write")
        write(991,*) "set term qt 1 noraise title 'Multigroup Diffusion Flux Profile'"
        if (mode == 5) then
            write(991,*) "set title 'Adjoint Flux (Importance Function) vs x for three energy groups (g_{3} < g_{1})'"
        else
            write(991,*) "set title 'Flux vs x for three energy groups'"
        end if
        write(991,*) "set xlabel 'x'"
        if (mode == 5) then
            write(991,*) "set ylabel 'Adjoint Flux (x)'"
        else
            write(991,*) "set ylabel 'phi(x)'"
        end if
        write(991,*) "plot 'flux.dat' using 1:2 with lines title 'g=1',\"
        write(991,*) "     'flux.dat' using 1:3 with lines title 'g=2',\"
        write(991,*) "     'flux.dat' using 1:4 with lines title 'g=3'"
        close(991)
        call system("env QT_QPA_PLATFORM=xcb /usr/bin/gnuplot -persist plot_flux.gp")

    end subroutine multigroup_diffusion_iter

end module m_driver


