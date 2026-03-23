program fem2d_main
    use omp_lib
    use m_constants
    use m_types
    use m_GMSH, only: ParseMesh, derive_case_nametag
    use m_quadrature, only: InitialiseQuadrature
    use m_finite_elements
    use m_sweep_order, only: InitiliaseConnectivityandNormals
    use m_material
    use m_outpVTK_dfem, only: export_dfem_vtk
    use m_transport, only: Transport_Sweep, Update_Scattering_Source_DGFEM, &
                           Update_Fission_Source_DGFEM, Calculate_Total_Production_DGFEM, &
                           Assemble_Fixed_Source_DGFEM, Precompute_Transport_Integrals
    implicit none

    real(dp), allocatable   :: angular_flux(:,:,:) ! (dof, n_angles, n_groups)
    real(dp), allocatable   :: scalar_flux(:,:), scalar_flux_old(:,:)
    real(dp), allocatable   :: scat_source(:,:), fiss_source(:,:), fixed_source(:,:)
    real(dp)                :: flux_diff

    type(t_mesh)            :: mesh
    type(t_quadrature)      :: Quad2D, Quad1D
    type(t_finite)          :: FE
    type(t_material), allocatable :: materials(:)
    type(t_sn_quadrature)   :: QuadSn 
    
    integer                 :: n_groups, n_dof
    integer                 :: outer_iter, max_outer_iter
    real(dp)                :: k_eff, k_eff_prime, total_prod_val
    logical                 :: is_eigenvalue_problem, is_adjoint
    character(len=32)       :: InputMesh

    InputMesh = "../input/trial_mesh.vtk"
    is_eigenvalue_problem  = .false.
    is_adjoint             = .false.
    n_groups               = 7
    QuadSn%order           = 16
    max_outer_iter         = 600

    call ParseMesh(InputMesh, FE, mesh, write_files = .true.)

    call InitialiseQuadrature(FE, mesh, Quad1D, Quad2D, QuadSn, is_adjoint)
    call InitialiseFiniteElements(FE, Quad1D, Quad2D)
    call InitialiseMaterials(materials, mesh, n_groups, "../MATS.txt", printout = .true.)
    call InitiliaseConnectivityandNormals(mesh, FE)

    call Precompute_Transport_Integrals(mesh, FE, Quad2D, Quad1D)

    n_dof = mesh%n_elems * FE%n_basis
    
    allocate(angular_flux(n_dof, QuadSn%NoAngles, n_groups))
    allocate(scalar_flux(n_dof, n_groups), scalar_flux_old(n_dof, n_groups))
    allocate(scat_source(n_dof, n_groups), fiss_source(n_dof, n_groups), fixed_source(n_dof, n_groups))

    angular_flux = 0.0_dp
    scalar_flux = 1.0_dp
    scalar_flux_old = 0.0_dp
    scat_source = 0.0_dp
    fiss_source = 0.0_dp
    fixed_source = 0.0_dp

    call Assemble_Fixed_Source_DGFEM(fixed_source, materials, mesh, FE, n_groups)

    write(*,*) merge(">>> Starting DGFEM Transport Power Iteration ... ", &
                     ">>> Starting DGFEM Transport Source Iteration ...", is_eigenvalue_problem)
    
    if (is_eigenvalue_problem) then
        call Calculate_Total_Production_DGFEM(total_prod_val, scalar_flux, materials, mesh, FE, n_groups, is_adjoint)
        if (total_prod_val > 1.0e-12_dp) then
            scalar_flux = scalar_flux / total_prod_val
        end if
    end if

    k_eff = 1.0_dp     
    do outer_iter = 1, max_outer_iter

        k_eff_prime = k_eff;            call Update_Fission_Source_DGFEM(fiss_source, scalar_flux, k_eff_prime, materials, mesh, FE, n_groups, is_adjoint)
        scalar_flux_old = scalar_flux;  call Update_Scattering_Source_DGFEM(scat_source, scalar_flux_old, materials, mesh, FE, n_groups, is_adjoint)
        scalar_flux = 0.0_dp;           call Transport_Sweep(mesh, FE, materials, QuadSn, angular_flux, scalar_flux, fixed_source, scat_source, fiss_source, n_groups)
                                    
        flux_diff = maxval(abs(scalar_flux - scalar_flux_old))
    
        if (is_eigenvalue_problem) then

            call Calculate_Total_Production_DGFEM(total_prod_val, scalar_flux, materials, mesh, FE, n_groups, is_adjoint)
            
            if (total_prod_val <= 1.0e-12_dp) then
                write(*,*) 'FATAL: Total production is zero or too small (no fission production).'
                write(*,*) '       Check material NuSigF values and mesh material assignments.'
                exit
            end if

            k_eff = total_prod_val * k_eff_prime

            scalar_flux = scalar_flux / total_prod_val
            
            write(*,'(A,I3,A,F12.8,A,ES12.4)') " Power Iteration: ", outer_iter, " k-eff: ", k_eff, "  Flux Diff: ", flux_diff
            if (abs(k_eff - k_eff_prime) < 1e-6) exit
        else
            write(*,'(A,I3,A,ES12.4)') " Source Iteration: ", outer_iter, " Flux Diff: ", flux_diff
            if (flux_diff < 1e-5) exit
        end if
    end do

    call export_dfem_vtk("../output/"//derive_case_nametag(InputMesh), mesh, FE, scalar_flux, n_groups)

    call system("gprof 2Diffusion > analysis.txt && python3 -m gprof2dot -f prof analysis.txt | dot -Tsvg -o ../output/profile_map"//derive_case_nametag(trim(InputMesh))//".svg") 
end program fem2d_main
