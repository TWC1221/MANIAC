program fem3d_main
    use omp_lib
    use m_constants
    use m_types
    use m_GMSH, only: ParseMesh, derive_case_nametag
    use m_quadrature, only: InitialiseQuadrature
    use m_finite_elements, only: InitialiseFiniteElements
    use m_sweep_order, only: InitialiseGeometry
    use m_material, only: InitialiseMaterials
    use m_utilities, only: print_splash, check_mesh_quality, print_timing_report, export_vtk!, PinPowerNormalisation
    use m_transport_precompute, only: InitialiseTransport
    use m_transport, only: Transport_Sweep, Calculate_Total_Production_DGFEM, Source_DGFEM
    implicit none

    ! --- Data Arrays ---
    real(dp), allocatable           :: angular_flux(:,:,:), scalar_flux(:,:), scalar_flux_old(:,:), total_src(:,:), pin_powers(:)

    ! ---   Objects   ---
    type(t_mesh)                    :: mesh
    type(t_quadrature)              :: Quad2D, Quad3D
    type(t_finite)                  :: FE
    type(t_material), allocatable   :: materials(:)
    type(t_sn_quadrature)           :: QuadSn 
    
    ! --- Sweep & Iteration Control ---
    integer, allocatable            :: sweep_order(:,:)
    integer                         :: n_groups, ref_ID(1)
    integer                         :: outer_iter, max_outer_iter, viewport_frequency, report_unit
    real(dp)                        :: k_eff, k_eff_prime, total_prod_val, flux_diff, tol
    logical                         :: is_eigenvalue_problem, is_adjoint, is_SEM, printout
    character(len=128)              :: InputMesh, report_path

    ! --- Configuration ---
    is_eigenvalue_problem  = .true.
    is_adjoint             = .false.
    is_SEM                 = .true.
    printout               = .true.
    FE%order               = 16
    QuadSn%order           = 16

    InputMesh              = "../input/trial_mesh.vtk"
    ref_ID                 = [103] !REFLECTIVE BOUNDARY ID
    n_groups               = 7

    max_outer_iter         = 2000
    viewport_frequency     = 75
    tol                    = 1e-6
    
    if (printout) call print_splash()
    
    call ParseMesh(InputMesh, FE, mesh, is_SEM)
    report_path = trim("../output/" // derive_case_nametag(InputMesh)) // "/performance_report" // trim(merge("_s","_f",is_SEM)) // &
    trim(merge("_adj","_fwd",is_adjoint)) // " n=" // trim(int_to_str(FE%order)) // " sn=" // trim(int_to_str(QuadSn%order)) // ".dat"
    open(newunit=report_unit, file=trim(report_path), status='replace')

    call InitialiseQuadrature(FE, Quad2D, Quad3D, QuadSn, is_adjoint, is_SEM) ! No printout, no file unit
    call InitialiseFiniteElements(FE, Quad2D, Quad3D) ! No printout, no file unit
    call InitialiseMaterials(materials, mesh, n_groups, "../input/MATS.txt", printout, file_unit=report_unit)
    call InitialiseGeometry(mesh, FE, QuadSn, sweep_order) ! No printout, no file unit
    call InitialiseTransport(mesh, FE, Quad3D, Quad2D, QuadSn, materials, n_groups) ! No printout, no file unit

    call check_mesh_quality(mesh, FE, printout, file_unit=report_unit)

    allocate(angular_flux(int(mesh%n_elems, 8) * FE%n_basis, QuadSn%n_angles, n_groups), &
             scalar_flux(int(mesh%n_elems, 8) * FE%n_basis, n_groups), &
             scalar_flux_old(int(mesh%n_elems, 8) * FE%n_basis, n_groups), &
             total_src(int(mesh%n_elems, 8) * FE%n_basis, n_groups))
    k_eff = 1.0_dp; angular_flux = 0.0; scalar_flux = 0.00001

    call timer_start('Power Iteration')
    do outer_iter = 1, max_outer_iter

        k_eff_prime = k_eff
        scalar_flux_old = scalar_flux

        call Source_DGFEM(total_src, scalar_flux, k_eff_prime, materials, mesh, FE, n_groups, is_adjoint, is_eigenvalue_problem)
        call Transport_Sweep(mesh, FE, QuadSn, angular_flux, scalar_flux, total_src, sweep_order, n_groups, ref_ID)
        
        flux_diff = maxval(abs(scalar_flux - scalar_flux_old))
        if (is_eigenvalue_problem) then
            call Calculate_Total_Production_DGFEM(total_prod_val, scalar_flux, materials, mesh, FE, n_groups, is_adjoint)
            k_eff = total_prod_val * k_eff_prime
            scalar_flux = scalar_flux / total_prod_val
            if (mod(outer_iter,10) == 0.0_dp) write(*,'(A,I3,A,F12.8,A,ES12.4)') " Outer: ", outer_iter, " Keff: ", k_eff
            if (abs(k_eff - k_eff_prime)/k_eff < tol) exit
        else
            write(*,'(A,I3,A,ES12.4)') " Outer: ", outer_iter, " Diff: ", flux_diff
            if (flux_diff < tol) exit
        end if

        if (mod(outer_iter,viewport_frequency) == 0.0_dp) call export_vtk("../output/"//derive_case_nametag(InputMesh), mesh, FE, QuadSN, scalar_flux, n_groups, pin_powers, is_SEM, is_adjoint, is_eigenvalue_problem, &
            max_outer_iter, outer_iter, tol, report_unit, k_eff, flux_diff)
    end do
    call timer_stop('Power Iteration')

    write(*,'(A,I3,A,ES12.4)') " >>> K_eff converged at iteration", outer_iter, " to", k_eff

    call export_vtk("../output/"//derive_case_nametag(InputMesh), mesh, FE, QuadSN, scalar_flux, n_groups, pin_powers, is_SEM, is_adjoint, is_eigenvalue_problem, &
        max_outer_iter, outer_iter, tol, report_unit, k_eff, flux_diff)
    call print_timing_report(printout, report_unit)

    close(report_unit)
    write(*,*) " >>> Output to "//derive_case_nametag(InputMesh)
end program fem3d_main