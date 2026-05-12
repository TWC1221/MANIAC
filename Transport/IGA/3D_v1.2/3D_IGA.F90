program fem3d_main
    use omp_lib
    use m_constants
    use m_types
    use m_asmg, only: read_asmg_mesh, derive_case_nametag, int_to_str
    use m_quadrature, only: InitialiseQuadrature
    use m_basis, only: InitialiseBasis
    use m_sweep_order, only: InitialiseGeometry
    use m_material, only: InitialiseMaterials
    use m_outpVTK_dfem, only: export_dfem_vtk
    use m_transport_precompute, only: InitialiseTransport, Precompute_Transport_Integrals
    use m_transport, only: Transport_Sweep, Calculate_Total_Production_DGFEM, Source_DGFEM
    implicit none

    ! --- Data Arrays ---
    real(dp), allocatable           :: angular_flux(:,:,:), scalar_flux(:,:), scalar_flux_old(:,:), total_src(:,:)

    ! ---   Objects   ---
    type(t_mesh), allocatable       :: mesh
    type(t_quadrature)              :: QuadFace, QuadVol
    type(t_finite)                  :: FE
    type(t_material), allocatable   :: materials(:)
    type(t_sn_quadrature)           :: QuadSn 
    
    ! --- Sweep Control ---
    integer, allocatable            :: sweep_order(:,:)
    integer                         :: ref_ID(1)
    integer                         :: outer_iter, max_outer_iter
    real(dp)                        :: k_eff, k_eff_prime, total_prod_val, flux_diff, tol, t_start, t_pre, t_solve, t_sweep_acc, t_src_acc, t1
    logical                         :: is_eigenvalue_problem, is_adjoint
    character(len=128)              :: InputMesh

    ! --- Configuration ---
    is_eigenvalue_problem  = .true.
    is_adjoint             = .false.

    InputMesh              = "../input/rod_test.asmg"
    ref_ID                 = [2] ! Match BC IDs used in .asmg files (e.g., reflective is 2)

    QuadSn%order           = 16

    max_outer_iter          = 900
    tol                    = 1e-7
    
    t_start = omp_get_wtime()

    call read_asmg_mesh(InputMesh, mesh) ! mesh is already allocated in the main program
    call InitialiseBasis(FE, mesh)
    call InitialiseQuadrature(FE, mesh, QuadFace, QuadVol, QuadSn, is_adjoint)
    call InitialiseMaterials(materials, mesh, "../input/MATS.txt", printout = .true.)
    call InitialiseGeometry(mesh, FE, QuadSn, QuadFace, QuadVol, sweep_order)
    call InitialiseTransport(mesh, FE, QuadSn, QuadVol, QuadFace, materials)

    allocate(angular_flux(mesh%n_elems * FE%n_basis, QuadSn%n_angles, mesh%n_groups), &
             scalar_flux(mesh%n_elems * FE%n_basis, mesh%n_groups), &
             scalar_flux_old(mesh%n_elems * FE%n_basis, mesh%n_groups), &
             total_src(mesh%n_elems * FE%n_basis, mesh%n_groups))

    t_pre = omp_get_wtime() - t_start
    k_eff = 1.0_dp; angular_flux = 0.0_dp; scalar_flux = 0.01_dp

    t_sweep_acc = 0.0_dp; t_src_acc = 0.0_dp
    t_solve = omp_get_wtime()

    do outer_iter = 1, max_outer_iter
        k_eff_prime = k_eff

        t1 = omp_get_wtime()
        call Source_DGFEM(total_src, scalar_flux, k_eff_prime, materials, mesh, FE, mesh%n_groups, is_adjoint, is_eigenvalue_problem)
        t_src_acc = t_src_acc + (omp_get_wtime() - t1)

        t1 = omp_get_wtime()
        call Transport_Sweep(mesh, FE, QuadSn, angular_flux, scalar_flux, total_src, sweep_order, ref_ID)
        t_sweep_acc = t_sweep_acc + (omp_get_wtime() - t1)
    
        if (is_eigenvalue_problem) then
            scalar_flux_old = scalar_flux
            call Calculate_Total_Production_DGFEM(total_prod_val, scalar_flux, materials, mesh, FE, is_adjoint)
            k_eff = total_prod_val * k_eff_prime
            
            if (total_prod_val > vsmall_number) then
                scalar_flux = scalar_flux / total_prod_val
                angular_flux = angular_flux / total_prod_val
            else
                write(*,*) "Warning: Fission production dropped to zero. Iteration halted to prevent overflow."
                exit
            end if

            ! Calculate difference between normalized fluxes for a valid convergence check
            flux_diff = maxval(abs(scalar_flux - scalar_flux_old))

            write(*,'(A,I6,A,F16.8,A,ES12.4)') " Outer: ", outer_iter, " k: ", k_eff, " Diff: ", flux_diff
            if (abs(k_eff - k_eff_prime)/k_eff < tol) exit
        else
            flux_diff = maxval(abs(scalar_flux - scalar_flux_old))
            scalar_flux_old = scalar_flux
            write(*,'(A,I3,A,ES12.4)') " Outer: ", outer_iter, " Diff: ", flux_diff
            if (flux_diff < tol) exit
        end if
    end do

    t_solve = omp_get_wtime() - t_solve

    write(*,'(/,A)') "--------------------------------------------------------"
    write(*,'(A)') "                 SOLVER TIMING SUMMARY                  "
    write(*,'(A)') "--------------------------------------------------------"
    write(*,'(A, F10.3, A)') "  Setup & Precompute: ", t_pre, " s"
    write(*,'(A, F10.3, A)') "  Iterative Solve:    ", t_solve, " s"
    write(*,'(A, F10.3, A)') "    - Total Sweeps:   ", t_sweep_acc, " s"
    write(*,'(A, F10.3, A)') "    - Total Sources:  ", t_src_acc, " s"
    write(*,'(A, F10.3, A)') "  Total Runtime:      ", omp_get_wtime() - t_start, " s"
    write(*,'(A)') "--------------------------------------------------------"

    call export_dfem_vtk("../output/"//derive_case_nametag(InputMesh), mesh, FE, QuadFace, QuadSN, &
                         scalar_flux, is_adjoint, refine_level_in=8, ang_flux_in=angular_flux)
    write(*,*) ">>> Simulation Complete."
end program