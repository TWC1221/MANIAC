program fem2d_main
    use omp_lib
    use m_constants
    use m_types
    use m_asmg, only: read_asmg_mesh
    use m_quadrature, only: InitialiseQuadrature
    use m_basis
    use m_sweep_order, only: InitialiseGeometry
    use m_material, only: InitialiseMaterials
    use m_constants, only: check_nan_scalar, check_nan_array, check_nan_matrix, int_to_str
    use m_outpVTK_dfem, only: export_dfem_vtk, PinPowerNormalisation, derive_case_nametag
    use m_transport_precompute, only: InitialiseTransport
    use m_transport, only: Transport_Sweep, Calculate_Total_Production_DGFEM, Source_DGFEM
    implicit none

    ! --- Data Arrays ---
    real(dp), allocatable           :: angular_flux(:,:,:), scalar_flux(:,:), scalar_flux_old(:,:), total_src(:,:), pin_powers(:)

    ! ---   Objects   ---
    type(t_mesh)                    :: mesh
    type(t_quadrature)              :: Quad2D, Quad1D
    type(t_finite)                  :: FE
    type(t_material), allocatable   :: materials(:)
    type(t_sn_quadrature)           :: QuadSn 
    
    ! --- Sweep Control ---
    integer, allocatable            :: sweep_order(:,:)
    integer                         :: n_groups, ref_ID(1)
    integer                         :: outer_iter, max_outer_iter
    real(dp)                        :: k_eff, k_eff_prime, total_prod_val, flux_diff, tol, t1, t2
    logical                         :: is_eigenvalue_problem, is_adjoint, is_SEM
    character(len=128)              :: InputMesh

    ! --- Configuration ---
    is_eigenvalue_problem  = .true.
    is_adjoint             = .false.

    InputMesh              = "../input/pincell_test.asmg"
    ref_ID                 = [-1] !B=101,R=102,T=103,L=104
    n_groups               = 7

    FE%order               = 2
    QuadSn%order           = 8

    max_outer_iter          = 600
    tol                    = 1e-4
    
    t1 = omp_get_wtime()

    call read_asmg_mesh(InputMesh, mesh)
    FE%n_basis = mesh%nloc

    call InitialiseQuadrature(FE, mesh, Quad1D, Quad2D, QuadSn, is_adjoint, is_SEM)
    call InitialiseFiniteElements(FE, Quad1D, Quad2D)
    call InitialiseMaterials(materials, mesh, n_groups, "../input/MATS.txt", printout = .true.)
    call InitialiseGeometry(mesh, FE, QuadSn, sweep_order)
    call InitialiseTransport(mesh, FE, Quad2D, Quad1D, QuadSn, materials, n_groups)

    allocate(angular_flux(mesh%n_elems * FE%n_basis, QuadSn%n_angles, n_groups), scalar_flux(mesh%n_elems * FE%n_basis, n_groups), scalar_flux_old(mesh%n_elems * FE%n_basis, n_groups), total_src(mesh%n_elems * FE%n_basis, n_groups))
    
    k_eff = 1.0_dp; angular_flux = 1.0_dp; scalar_flux = 1.0_dp
    call check_nan_scalar(k_eff, "k_eff (initial)", "main program")
    call check_nan_matrix(scalar_flux, "scalar_flux (initial)", "main program")

    do outer_iter = 1, max_outer_iter
        k_eff_prime = k_eff
        scalar_flux_old = scalar_flux
        print*,sum(scalar_flux), "Total Flux at Outer Iteration ", outer_iter
        call Source_DGFEM(total_src, scalar_flux, k_eff_prime, materials, mesh, FE, n_groups, is_adjoint, is_eigenvalue_problem)
        print*,sum(total_src), "Total Source at Outer Iteration ", outer_iter
        call Transport_Sweep(mesh, FE, QuadSn, angular_flux, scalar_flux, total_src, sweep_order, n_groups, ref_ID)
    
        flux_diff = maxval(abs(scalar_flux - scalar_flux_old))
    
        call check_nan_scalar(flux_diff, "flux_diff", "main program, Outer Iter "//int_to_str(outer_iter))
        if (is_eigenvalue_problem) then
            call Calculate_Total_Production_DGFEM(total_prod_val, scalar_flux, materials, mesh, FE, n_groups, is_adjoint)
            
            if (isnan(total_prod_val)) then
                write(*,*) "FATAL: NaN detected in total_prod_val at outer iteration ", outer_iter, ". Simulation stopped."
                stop "NaN in total production."
            else if (total_prod_val < dp_EPSILON) then ! Using dp_EPSILON from m_constants for a more robust check
                write(*,*) "FATAL: Total production (", total_prod_val, ") is too small (< ", dp_EPSILON, ") at outer iteration ", outer_iter, ". Current k_eff_prime = ", k_eff_prime, ". Simulation stopped."
                stop "Total production too small."
            end if

            ! The check_nan_scalar for total_prod_val is now redundant as isnan is explicitly handled above.
            ! call check_nan_scalar(total_prod_val, "total_prod_val", "main program, Outer Iter "//int_to_str(outer_iter))
            k_eff = total_prod_val * k_eff_prime

            if (isnan(k_eff) .or. k_eff < dp_EPSILON) then
                write(*,*) "FATAL: k_eff (", k_eff, ") became NaN or too small (< ", dp_EPSILON, ") at outer iteration ", outer_iter, ". Simulation stopped."
                stop
            end if
            scalar_flux = scalar_flux / total_prod_val
            write(*,'(A,I3,A,F12.8,A,ES12.4)') " Outer: ", outer_iter, " k: ", k_eff, " Diff: ", flux_diff
            if (abs(k_eff - k_eff_prime)/k_eff < tol) exit
        else
            write(*,'(A,I3,A,ES12.4)') " Outer: ", outer_iter, " Diff: ", flux_diff
            if (flux_diff < tol) exit
        end if
    end do

   ! if (index(derive_case_nametag(InputMesh),"C5G7") > 0) call PinPowerNormalisation(mesh, FE, materials, Quad2D, n_groups, scalar_flux, angular_flux, total_src, pin_powers, "../input/PINS.txt")
    
    call export_dfem_vtk("../output/"//derive_case_nametag(InputMesh), mesh, FE, QuadSN, scalar_flux, n_groups, pin_powers, is_SEM, is_adjoint)
    write(*,*) ">>> Simulation Complete."

    t2 = omp_get_wtime()

    print*,t2-t1

end program fem2d_main