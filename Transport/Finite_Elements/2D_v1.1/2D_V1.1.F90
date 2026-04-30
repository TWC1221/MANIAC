program fem2d_main
    use omp_lib
    use m_constants
    use m_types
    use m_GMSH, only: ParseMesh, derive_case_nametag
    use m_quadrature, only: InitialiseQuadrature
    use m_finite_elements, only: InitialiseFiniteElements
    use m_sweep_order, only: InitialiseGeometry
    use m_material, only: InitialiseMaterials
    use m_outpVTK_dfem, only: export_dfem_vtk, PinPowerNormalisation
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
    integer                         :: n_groups, ref_ID(2)
    integer                         :: outer_iter, max_outer_iter
    real(dp)                        :: k_eff, k_eff_prime, total_prod_val, flux_diff, tol, t1, t2
    logical                         :: is_eigenvalue_problem, is_adjoint, is_SEM
    character(len=128)              :: InputMesh

    ! --- Configuration ---
    is_eigenvalue_problem  = .true.
    is_adjoint             = .false.
    is_SEM                 = .true.

    InputMesh              = "../input/C5G7vol_1.vtk"
    ref_ID                 = [103, 104] !B=101,R=102,T=103,L=104
    n_groups               = 7

    FE%order               = 1
    QuadSn%order           = 8

    max_outer_iter          = 600
    tol                    = 1e-6
    
    t1 = omp_get_wtime()

    call ParseMesh(InputMesh, FE, mesh, is_SEM, write_files = .true.)

    call InitialiseQuadrature(FE, mesh, Quad1D, Quad2D, QuadSn, is_adjoint, is_SEM)
    call InitialiseFiniteElements(FE, Quad1D, Quad2D)
    call InitialiseMaterials(materials, mesh, n_groups, "../input/MATS.txt", printout = .true.)
    call InitialiseGeometry(mesh, FE, QuadSn, sweep_order)
    print*,sweep_order(1,1:8)
    call InitialiseTransport(mesh, FE, Quad2D, Quad1D, QuadSn, materials, n_groups)

    allocate(angular_flux(mesh%n_elems * FE%n_basis, QuadSn%n_angles, n_groups), scalar_flux(mesh%n_elems * FE%n_basis, n_groups), scalar_flux_old(mesh%n_elems * FE%n_basis, n_groups), total_src(mesh%n_elems * FE%n_basis, n_groups))
    k_eff = 1.0_dp; angular_flux = 0.0; scalar_flux = 0.001

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
            write(*,'(A,I3,A,F16.8,A,ES12.4)') " Outer: ", outer_iter, " k: ", k_eff, " Diff: ", flux_diff
            if (abs(k_eff - k_eff_prime)/k_eff < tol) exit
        else
            write(*,'(A,I3,A,ES12.4)') " Outer: ", outer_iter, " Diff: ", flux_diff
            if (flux_diff < tol) exit
        end if
    end do

    if (index(derive_case_nametag(InputMesh),"C5G7") > 0) call PinPowerNormalisation(mesh, FE, materials, Quad2D, n_groups, scalar_flux, angular_flux, total_src, pin_powers, "../input/PINS.txt")
    
    call export_dfem_vtk("../output/"//derive_case_nametag(InputMesh), mesh, FE, QuadSN, scalar_flux, n_groups, pin_powers, is_SEM, is_adjoint)
    write(*,*) ">>> Simulation Complete."

    t2 = omp_get_wtime()

    print*,t2-t1

end program fem2d_main