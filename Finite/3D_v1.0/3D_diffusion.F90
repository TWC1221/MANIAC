#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

program fem3d_main
    USE OMP_LIB
    use m_constants
    use m_types
    use m_GMSH
    use m_quadrature
    use m_finite_elements
    use m_material
    use m_boundaries
    use m_outpVTK
    use m_PCG
    use m_petsc
    use m_timing

    implicit none

    Mat, allocatable        :: A_MAT(:)    
    Mat, allocatable        :: MAT_F_PETSC(:,:), MAT_S_PETSC(:,:)
    Vec, allocatable        :: B_VEC(:), X_VEC(:), X_PRIME_VEC(:)
    Vec, allocatable        :: FixedSrc_PETSC(:), PROD_VEC(:), temp_vec(:)
    KSP, allocatable        :: ksp_solvers(:)
    PetscErrorCode          :: ierr

    type(t_mat), allocatable :: A_PCG(:)
    type(t_mat), allocatable :: MAT_F_PCG(:,:), MAT_S_PCG(:,:)
    type(t_vec), allocatable :: B_PCG(:), X_PCG(:), X_PRIME_PCG(:), FixedSrc_PCG(:)

    type(t_mesh)            :: mesh
    type(t_bc_config)       :: bc_config(6)
    type(t_quadrature)      :: Quad, QuadBound
    type(t_finite)          :: FE
    type(t_material), allocatable :: materials(:)
    
    integer                 :: n_groups, i, k, outer_iter, omp_num_threads
    real(dp)                :: max_phi_change, group_diff
    integer                 :: solver_choice, preconditioner_choice
    integer                 :: max_outer_iter, max_CG_iter
    logical                 :: is_eigenvalue_problem, is_adjoint
    real(dp)                :: k_eff, k_eff_prime, total_production
    character(len=32)       :: InputMesh
    logical                 :: verbose

    InputMesh = "../input/C5G7v9.vtk"
    is_eigenvalue_problem   = .true.
    is_adjoint              = .false.
    n_groups               = 7

    solver_choice          = SOLVER_KSP_CG
    preconditioner_choice  = PRECON_NONE
    omp_num_threads        = 40
    max_outer_iter         = 500          
    max_CG_iter            = 100000          
    verbose                = .true.

    call timer_start('Total Execution')
    call timer_start('Initialization')

    call timer_start('PETSc Init')
    if (solver_choice /= SOLVER_PCG) then
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        call omp_set_num_threads(omp_num_threads)
    end if
    call timer_stop('PETSc Init')
    
    call timer_start('Mesh Parsing')
    call parse_GMSH(InputMesh, mesh, write_files = .true., verbose=.true.)
    call timer_stop('Mesh Parsing')

    FE%n_basis = mesh%nloc
    FE%order   = nint(real(FE%n_basis, dp)**(1.0_dp/3.0_dp)-1.0_dp)

    call timer_start('FE Generation')
    call QuadrilateralQuadrature(QuadBound, FE%order + 1)
    call HexahedralQuadrature(Quad, FE%order + 1, writeout = .false.) 

    call InitialiseFiniteElements(FE, Quad, QuadBound, write_nodal_map = (.false.))
    call timer_stop('FE Generation')

    call timer_start('Material Init')
    call InitialiseMaterials(materials, mesh, n_groups, "../MATS.txt", printout = (verbose))
    call timer_stop('Material Init')

    call timer_start('Boundary Init')
    call InitialiseBoundaries(bc_config(1), 101, BC_REFLECTIVE, 0.0_dp) !XMIN
    call InitialiseBoundaries(bc_config(2), 102, BC_VACUUM, 0.0_dp) !XMAX
    call InitialiseBoundaries(bc_config(3), 103, BC_REFLECTIVE, 0.0_dp) !YMIN
    call InitialiseBoundaries(bc_config(4), 104, BC_VACUUM, 0.0_dp) !YMAX
    call InitialiseBoundaries(bc_config(5), 105, BC_REFLECTIVE, 0.0_dp) !ZMIN
    call InitialiseBoundaries(bc_config(6), 106, BC_VACUUM, 0.0_dp) !ZMAX
    call timer_stop('Boundary Init')

    call timer_stop('Initialization')

    call timer_start('Matrix Assembly')
    if (solver_choice /= SOLVER_PCG) then
        if (verbose) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-mat_view", "draw", ierr)
        ! if (verbose) call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-draw_pause", "-1", ierr)
        allocate(A_MAT(n_groups), B_VEC(n_groups), X_VEC(n_groups), X_PRIME_VEC(n_groups), PROD_VEC(n_groups), FixedSrc_PETSC(n_groups), temp_vec(n_groups), ksp_solvers(n_groups))
        do i = 1, n_groups
            call VecCreate(PETSC_COMM_SELF, B_VEC(i), ierr)
            call VecSetSizes(B_VEC(i), PETSC_DECIDE, mesh%n_nodes, ierr)
            call VecSetFromOptions(B_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), X_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), X_PRIME_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), temp_vec(i), ierr)
            call VecSet(X_VEC(i), merge(0.001_dp, 0.0_dp, is_eigenvalue_problem), ierr)
            call VecSet(X_PRIME_VEC(i), merge(0.001_dp, 0.0_dp, is_eigenvalue_problem), ierr)  
            call KSPCreate(PETSC_COMM_SELF, ksp_solvers(i), ierr)
        end do
        if (verbose) write(*,'(A)') ">>> Assembling Global Matrix..."
        call ASSEMBLEMultigroupMAT_PETSc(A_MAT, MAT_F_PETSC, MAT_S_PETSC, PROD_VEC, FixedSrc_PETSC, mesh, FE, Quad, materials, n_groups, is_adjoint, printout=verbose)
    else 
        allocate(A_PCG(n_groups), B_PCG(n_groups), X_PCG(n_groups), X_PRIME_PCG(n_groups))
        do i = 1, n_groups
            call PCG_VEC_INIT(mesh%n_nodes, X_PCG(i))
            call PCG_VEC_INIT(mesh%n_nodes, X_PRIME_PCG(i))
            call PCG_VEC_INIT(mesh%n_nodes, B_PCG(i))
            X_PCG(i)%vec = merge(0.001_dp, 0.0_dp, is_eigenvalue_problem)
            X_PRIME_PCG(i)%vec = merge(0.001_dp, 0.0_dp, is_eigenvalue_problem)
        end do
        if (verbose) write(*,'(A)') ">>> Assembling Global Matrix..."
        call assemble_PCG_matrix(materials, A_PCG, mesh, FE, Quad, n_groups, mesh%mats)
        call assemble_source_matrices_pcg(MAT_F_PCG, MAT_S_PCG, FixedSrc_PCG, mesh, FE, Quad, materials, n_groups, is_adjoint)
    end if
    call timer_stop('Matrix Assembly')

    call timer_start('KSP + Boundary Setup')
    do i = 1, n_groups
        do k = 1, size(bc_config)
            if (solver_choice /= SOLVER_PCG) then
                call apply_bcs(mesh, FE, Quadbound, bc_config(k), A_petsc=A_MAT(i), printout = .false.)
            else
                call apply_bcs(mesh, FE, Quadbound, bc_config(k), A_pcg=A_PCG(i), printout = .false.)
            end if
        end do

        if (solver_choice /= SOLVER_PCG) then
            call SetupKSP_PETSc(ksp_solvers(i), A_MAT(i), solver_choice, preconditioner_choice)
        end if
    end do
    call timer_stop('KSP + Boundary Setup')

    if (verbose .and. is_eigenvalue_problem) write(*,'(A)') ">>> Starting Power Iteration..."
    if (verbose .and. (.not. is_eigenvalue_problem)) write(*,'(A)') ">>> Solving for Neutron Flux..."

    call timer_start('Solver Loop')

    k_eff = 1.0_dp     
    do outer_iter = 1, max_outer_iter
        k_eff_prime = k_eff
        max_phi_change = 0.0_dp

        call timer_start('Source Evaluation')
        do i = 1, n_groups
            if (solver_choice /= SOLVER_PCG) then
                call EVALMultigroup_PETSc(B_VEC(i), MAT_F_PETSC, MAT_S_PETSC, FixedSrc_PETSC, &
                                                      X_VEC, X_PRIME_VEC, n_groups, k_eff_prime, i, temp_vec(i), is_eigenvalue_problem)
            else
                call assemble_multigroup_source_pcg(B_PCG(i), MAT_F_PCG, MAT_S_PCG, FixedSrc_PCG, &
                                                    X_PCG, X_PRIME_PCG, n_groups, k_eff_prime, i, is_eigenvalue_problem)
            end if
        end do
        call timer_stop('Source Evaluation')

        call timer_start('Linear Solve')
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, group_diff, ierr) REDUCTION(max:max_phi_change)
        do i = 1, n_groups
            if (solver_choice /= SOLVER_PCG) then
                call KSPSolve(ksp_solvers(i), B_VEC(i), X_VEC(i), ierr)
                call VecWAXPY(B_VEC(i), -1.0_dp, X_PRIME_VEC(i), X_VEC(i), ierr)
                call VecNorm(B_VEC(i), NORM_INFINITY, group_diff, ierr)
                max_phi_change = max(max_phi_change, group_diff)
            else
                call CGSolve(A_PCG(i), X_PCG(i), B_PCG(i), preconditioner_choice, max_CG_iter)
                max_phi_change = max(max_phi_change, maxval(abs(X_PCG(i)%vec - X_PRIME_PCG(i)%vec)))
            end if
        end do
        !$OMP END PARALLEL DO
        call timer_stop('Linear Solve')

        call timer_start('Eigen Update')
        if (is_eigenvalue_problem) then
            if (solver_choice /= SOLVER_PCG) then
                call EVALProd_PETSc(total_production, X_VEC, PROD_VEC, n_groups)
                do i = 1, n_groups
                    call VecScale(X_VEC(i), 1.0_dp / total_production, ierr)
                    call VecCopy(X_VEC(i), X_PRIME_VEC(i), ierr)
                end do
            else
                call calculate_total_production_pcg(total_production, X_PCG, mesh, FE, Quad, materials, n_groups)
                do i = 1, n_groups
                    X_PCG(i)%vec = X_PCG(i)%vec / total_production
                    X_PRIME_PCG(i)%vec = X_PCG(i)%vec
                end do
            end if

            k_eff = k_eff_prime * total_production

            if (verbose) write(*,'(A,I3,A,F12.8,A,ES12.4)') " [ SOLVER ] :: Iteration: ", outer_iter, " k-eff: ", k_eff, "  dPhi: ", max_phi_change
            if (abs(k_eff - k_eff_prime) < 1e-5 .and. max_phi_change < 1e-5) exit
        else
            if (solver_choice /= SOLVER_PCG) then
                do i = 1, n_groups
                    call VecCopy(X_VEC(i), X_PRIME_VEC(i), ierr)
                end do
            else
                do i = 1, n_groups
                    X_PRIME_PCG(i)%vec = X_PCG(i)%vec
                end do
            end if
            if (verbose) write(*,'(A,I3,A,ES12.4)') " [ SOLVER ] :: Source Iteration: ", outer_iter, " Max Change: ", max_phi_change
            if (max_phi_change < 1e-8) exit
        end if 
        call timer_stop('Eigen Update')
    end do

    call timer_stop('Solver Loop')

    call timer_start('Output')

    if (solver_choice /= SOLVER_PCG) then
        call export_vtk_petsc("../output/"//derive_case_nametag(InputMesh), FE, mesh, X_VEC, n_groups, 2, is_adjoint, printout=verbose)
        call PetscFinalize(ierr)
    else
        call export_vtk_pcg("../output/"//derive_case_nametag(InputMesh), FE, mesh, X_PCG, n_groups, 2, .false., printout=verbose)
    end if

    call timer_stop('Output')

    call timer_stop('Total Execution')
    
    call print_timing_report(verbose)
    call system("gprof 3Diffusion > analysis.txt && python3 -m gprof2dot -f prof analysis.txt | dot -Tsvg -o ../output/profile_map"//derive_case_nametag(trim(InputMesh))//".svg") 

end program fem3d_main