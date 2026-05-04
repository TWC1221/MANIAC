program fem2d_main
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
    type(t_mesh)                    :: mesh
    type(t_quadrature)              :: Quad2D, Quad1D
    type(t_finite)                  :: FE
    type(t_material), allocatable   :: materials(:)
    type(t_sn_quadrature)           :: QuadSn 
    
    ! --- Sweep Control ---
    integer, allocatable            :: sweep_order(:,:)
    integer                         :: ref_ID(1)
    integer                         :: outer_iter, max_outer_iter
    real(dp)                        :: k_eff, k_eff_prime, total_prod_val, flux_diff, tol, t1, t2
    logical                         :: is_eigenvalue_problem, is_adjoint
    character(len=128)              :: InputMesh

    ! --- Configuration ---
    is_eigenvalue_problem  = .true.
    is_adjoint             = .false.

    InputMesh              = "../input/quatre.asmg" !squares_IGA_FDG
    ref_ID                 = [2]

    QuadSn%order           =  16

    max_outer_iter          = 900
    tol                    = 1e-7
    
    t1 = omp_get_wtime()

    call read_asmg_mesh(InputMesh, mesh)
    call InitialiseBasis(FE, mesh)
    call InitialiseQuadrature(FE, mesh, Quad1D, Quad2D, QuadSn, is_adjoint)
    call InitialiseMaterials(materials, mesh, "../input/MATS.txt", printout = .true.)
    call InitialiseGeometry(mesh, FE, QuadSn, Quad2D, Quad1D, sweep_order)
    call InitialiseTransport(mesh, FE, QuadSn, Quad2D, Quad1D, materials)

    allocate(angular_flux(mesh%n_elems * FE%n_basis, QuadSn%n_angles, mesh%n_groups), &
             scalar_flux(mesh%n_elems * FE%n_basis, mesh%n_groups), &
             scalar_flux_old(mesh%n_elems * FE%n_basis, mesh%n_groups), &
             total_src(mesh%n_elems * FE%n_basis, mesh%n_groups))
    k_eff = 1.0_dp; angular_flux = 0.0; scalar_flux = 0.01_dp

    do outer_iter = 1, max_outer_iter

        k_eff_prime = k_eff
        scalar_flux_old = scalar_flux

        call Source_DGFEM(total_src, scalar_flux, k_eff_prime, materials, mesh, FE, mesh%n_groups, is_adjoint, is_eigenvalue_problem)
        call Transport_Sweep(mesh, FE, QuadSn, angular_flux, scalar_flux, total_src, sweep_order, ref_ID)
    
        flux_diff = maxval(abs(scalar_flux - scalar_flux_old))
    
        if (is_eigenvalue_problem) then
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
            write(*,'(A,I3,A,ES12.4)') " Outer: ", outer_iter, " Diff: ", flux_diff
            if (flux_diff < tol) exit
        end if
    end do

    call export_dfem_vtk("../output/"//derive_case_nametag(InputMesh), mesh, FE, Quad1D, QuadSN, scalar_flux, is_adjoint, refine_level_in=8)
    write(*,*) ">>> Simulation Complete."

    t2 = omp_get_wtime()

    print*,t2-t1

end program fem2d_main