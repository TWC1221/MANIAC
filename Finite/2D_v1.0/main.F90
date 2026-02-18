#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

program fem2d_main
    use m_constants
    use m_types
    use m_GMSH,            only: parse
    use m_quadrature
    use m_finite_elements
    use m_material
    use m_boundaries
    use m_outpVTK
    use m_PCG
    use m_petsc

    implicit none

    Mat, allocatable        :: A_MAT(:)    
    Vec, allocatable        :: B_VEC(:), X_VEC(:)    
    KSP, allocatable        :: ksp_solvers(:)
    PetscErrorCode          :: ierr

    type(t_mat), allocatable :: A_PCG(:)
    type(t_vec), allocatable :: B_PCG(:), X_PCG(:)

    type(t_mesh)            :: mesh
    type(t_bc_config)       :: bc_config(2)
    type(t_quadrature)      :: Quad, QuadBound
    type(t_finite)          :: FE
    type(t_material), allocatable :: materials(:)
    
    integer                 :: n_groups, i, k, outer_iter
    real(dp), allocatable   :: S_ext(:,:)
    integer                 :: solver_choice, max_iter_inner, max_iter_power
    logical                 :: is_eigenvalue_problem
    real(dp)                :: k_eff, k_old, norm_phi, total_prod
    integer, parameter      :: SOLVER_PCG = 1, SOLVER_PETSC = 2

    solver_choice          = SOLVER_PCG
    is_eigenvalue_problem  = .true.         
    n_groups               = 1
    max_iter_power         = 100           
    max_iter_inner         = 10000          
    k_eff                  = 1.0_dp        
    
    write(*,*) ">>> Initialising Problem ... ", merge("EIGENVALUE", "FIXED     ",is_eigenvalue_problem) 

    call parse("../VTKmesh.vtk", mesh, .false.)
    FE%n_basis = mesh%nloc
    FE%order   = nint(sqrt(real(FE%n_basis, dp))) - 1

    call Get1DLineQuad(FE%order + 1, QuadBound)
    call QuadrilateralQuadrature(Quad, FE%order + 1)

    call InitialiseFE(FE, Quad, QuadBound)
    call InitialiseMaterials(materials, n_groups)

    call InitialiseBoundaries(bc_config(1), 2, BC_VACUUM, 0.0_dp)
    call InitialiseBoundaries(bc_config(2), 4, BC_REFLECTIVE, 0.0_dp)

    if (solver_choice == SOLVER_PETSC) then
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        allocate(A_MAT(n_groups), B_VEC(n_groups), X_VEC(n_groups), ksp_solvers(n_groups))
    else
        allocate(A_PCG(n_groups), B_PCG(n_groups), X_PCG(n_groups))
    end if

    allocate(S_ext(mesh%n_elems, n_groups))
    S_ext = 0.0_dp
    if (.not. is_eigenvalue_problem) then
        do i = 1, mesh%n_elems
            if (mesh%mats(i) == 1) S_ext(i, :) = 10.0_dp
        end do
    end if

    if (solver_choice == SOLVER_PETSC) then
        do i = 1, n_groups
            call VecCreate(PETSC_COMM_WORLD, B_VEC(i), ierr)
            call VecSetSizes(B_VEC(i), PETSC_DECIDE, mesh%n_nodes, ierr)
            call VecSetFromOptions(B_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), X_VEC(i), ierr)
            call VecSet(X_VEC(i), 1.0_dp, ierr) 
            call KSPCreate(PETSC_COMM_WORLD, ksp_solvers(i), ierr)
        end do
    else
        do i = 1, n_groups
            call PCG_VEC_INIT(mesh%n_nodes, X_PCG(i))
            call PCG_VEC_INIT(mesh%n_nodes, B_PCG(i))
            X_PCG(i)%vec = 1.0_dp
        end do
    end if

    !=================================================================================================================

    write(*,*) merge(">>> Starting Power Iteration ...", &
                     ">>> Starting Global Assembly ...", is_eigenvalue_problem)

    do outer_iter = 1, merge(max_iter_power, 1, is_eigenvalue_problem)
        k_old = k_eff
        
        do i = 1, n_groups
            if (solver_choice == SOLVER_PETSC) then
                call assemble_petsc_matrix(materials, A_MAT, mesh, FE, Quad, n_groups, mesh%mats)
            else
                call assemble_PCG_matrix(materials, A_PCG, mesh, FE, Quad, n_groups, mesh%mats)
            end if

            if (is_eigenvalue_problem) then
                if (solver_choice == SOLVER_PETSC) then
                    call assemble_multigroup_fission_vec(B_VEC(i), mesh, FE, Quad, materials, X_VEC, n_groups, k_eff)
                else
                    call assemble_multigroup_fission_pcg(B_PCG(i), mesh, FE, Quad, materials, X_PCG, n_groups, k_eff)
                end if
            else
                if (solver_choice == SOLVER_PETSC) then
                    call assemble_petsc_source_vec(B_VEC, mesh, FE, Quad, n_groups, S_ext)
                else
                    call assemble_PCG_source_vec(B_PCG, mesh, FE, Quad, n_groups, S_ext)
                end if
            end if

            do k = 1, size(bc_config)
                if (solver_choice == SOLVER_PETSC) then
                    call apply_bcs(mesh, bc_config(k), A_petsc=A_MAT(i), B_petsc=B_VEC(i), &
                                   solver_choice=SOLVER_PETSC)
                    call MatAssemblyBegin(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
                    call MatAssemblyEnd(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
                    call KSPSetOperators(ksp_solvers(i), A_MAT(i), A_MAT(i), ierr)
                else
                    call apply_bcs(mesh, bc_config(k), A_pcg=A_PCG(i), B_pcg=B_PCG(i), &
                                   solver_choice=SOLVER_PCG)
                end if
            end do

            if (solver_choice == SOLVER_PETSC) then 
                call KSPSolve(ksp_solvers(i), B_VEC(i), X_VEC(i), ierr)
            else
                call CGSolve(A_PCG(i), X_PCG(i), B_PCG(i), PCG_PRECON_NONE, max_iter_inner)
            end if
        end do

        if (.not. is_eigenvalue_problem) exit 

        if (solver_choice == SOLVER_PETSC) then
            call calculate_total_production_petsc(total_prod, X_VEC, mesh, FE, Quad, materials, n_groups)
        else
            call calculate_total_production_pcg(total_prod, X_PCG, mesh, FE, Quad, materials, n_groups)
        end if

        k_eff = total_prod 

        write(*,'(A,I3,A,F12.8)') " Iteration: ", outer_iter, "  k-eff: ", k_eff
        if (abs(k_eff - k_old) < 1e-7) exit

        do i = 1, n_groups
            if (solver_choice == SOLVER_PETSC) then
                call VecScale(X_VEC(i), 1.0_dp / k_eff, ierr)
            else
                X_PCG(i)%vec = X_PCG(i)%vec / k_eff
            end if
        end do

    end do

    if (solver_choice == SOLVER_PETSC) then
        call export_vtk_petsc("solution_output", FE, mesh, X_VEC, n_groups, 2, .true.)
        do i = 1, n_groups
            call KSPDestroy(ksp_solvers(i), ierr)
            call MatDestroy(A_MAT(i), ierr)
            call VecDestroy(B_VEC(i), ierr)
            call VecDestroy(X_VEC(i), ierr)
        end do
        call PetscFinalize(ierr)
    else
        call export_vtk_pcg("solution_output", FE, mesh, X_PCG, n_groups, 2, .true.)
    end if
    if (allocated(S_ext)) deallocate(S_ext)
    write(*,*) ">>> Simulation Finished."

    call system("gprof exe > analysis.txt && python3 -m gprof2dot -f prof analysis.txt | dot -Tsvg -o profile_map.svg") 

end program fem2d_main