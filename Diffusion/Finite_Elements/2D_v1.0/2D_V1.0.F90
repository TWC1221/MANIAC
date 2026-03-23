#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

program fem2d_main
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

    implicit none

    Mat, allocatable        :: A_MAT(:)    
    Vec, allocatable        :: B_VEC(:), X_VEC(:), X_PRIME_VEC(:)    
    KSP, allocatable        :: ksp_solvers(:)
    PetscErrorCode          :: ierr

    type(t_mat), allocatable :: A_PCG(:)
    type(t_vec), allocatable :: B_PCG(:), X_PCG(:), X_PRIME_PCG(:)

    type(t_mesh)            :: mesh
    type(t_bc_config)       :: bc_config(4)
    type(t_quadrature)      :: Quad, QuadBound
    type(t_finite)          :: FE
    type(t_material), allocatable :: materials(:)
    
    integer                 :: n_groups, i, k, inner_iter, outer_iter
    real(dp), allocatable   :: S_ext(:,:)
    integer                 :: solver_choice, preconditioner_choice
    integer                 :: max_inner_iter, max_outer_iter, max_CG_iter
    logical                 :: is_eigenvalue_problem, is_adjoint
    real(dp)                :: k_eff, k_eff_prime
    integer, parameter      :: SOLVER_PCG = 1, SOLVER_PETSC = 2

    is_eigenvalue_problem   = .true.
    is_adjoint              = .false.
    solver_choice          = SOLVER_PCG
    preconditioner_choice  = PRECON_NONE
    n_groups               = 7
    max_outer_iter         = 500          
    max_inner_iter         = 1         
    max_CG_iter            = 100000          

    call parse_GMSH("../input/C5G7surf_1_1_f_p.vtk", mesh, write_files = .true.)
    FE%n_basis = mesh%nloc
    FE%order   = nint(sqrt(real(FE%n_basis, dp))) - 1

    call Get1DLineQuad(FE%order + 1, QuadBound)
    call QuadrilateralQuadrature(Quad, FE%order + 1)

    call InitialiseFiniteElements(FE, Quad, QuadBound, write_nodal_map = .false.)
    call InitialiseMaterials(materials, mesh, n_groups, "../MATS.txt", printout = .true.)

    call InitialiseBoundaries(bc_config(1), 101, BC_DIRICHLET, 0.0_dp)
    call InitialiseBoundaries(bc_config(2), 101, BC_DIRICHLET, 0.0_dp) 
    call InitialiseBoundaries(bc_config(3), 102, BC_REFLECTIVE, 0.0_dp)
    call InitialiseBoundaries(bc_config(4), 102, BC_REFLECTIVE, 0.0_dp) 

    allocate(S_ext(mesh%n_elems, n_groups)); S_ext = 1.0_dp

    write(*,*) ">>> Allocating Memory ..."
    if (solver_choice == SOLVER_PETSC) then
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        allocate(A_MAT(n_groups), B_VEC(n_groups), X_VEC(n_groups), X_PRIME_VEC(n_groups), ksp_solvers(n_groups))
        do i = 1, n_groups
            call VecCreate(PETSC_COMM_WORLD, B_VEC(i), ierr)
            call VecSetSizes(B_VEC(i), PETSC_DECIDE, mesh%n_nodes, ierr)
            call VecSetFromOptions(B_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), X_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), X_PRIME_VEC(i), ierr)
            call VecSet(X_VEC(i), merge(1.0_dp, 0.0_dp, is_eigenvalue_problem), ierr)
            call VecSet(X_PRIME_VEC(i), merge(1.0_dp, 0.0_dp, is_eigenvalue_problem), ierr)  
            call KSPCreate(PETSC_COMM_WORLD, ksp_solvers(i), ierr)
        end do
        call assemble_petsc_matrix(materials, A_MAT, mesh, FE, Quad, n_groups, mesh%mats)
    else
        allocate(A_PCG(n_groups), B_PCG(n_groups), X_PCG(n_groups), X_PRIME_PCG(n_groups))
        do i = 1, n_groups
            call PCG_VEC_INIT(mesh%n_nodes, X_PCG(i))
            call PCG_VEC_INIT(mesh%n_nodes, B_PCG(i))
            call PCG_VEC_INIT(mesh%n_nodes, X_PRIME_PCG(i))
            X_PCG(i)%vec = merge(1.0_dp, 0.0_dp, is_eigenvalue_problem)
            X_PRIME_PCG(i)%vec = merge(1.0_dp, 0.0_dp, is_eigenvalue_problem)
        end do
        call assemble_PCG_matrix(materials, A_PCG, mesh, FE, Quad, n_groups, mesh%mats)
    end if

    write(*,*) ">>> Assembling Global Matrix ..."
    do i = 1, n_groups
        do k = 1, size(bc_config)
            if (solver_choice == SOLVER_PETSC) then
                call apply_bcs(mesh, FE, Quadbound, bc_config(k), A_petsc=A_MAT(i), solver_choice=SOLVER_PETSC)
            else
                call apply_bcs(mesh, FE, Quadbound, bc_config(k), A_pcg=A_PCG(i), solver_choice=SOLVER_PCG)
            end if
        end do
        if (solver_choice == SOLVER_PETSC) then
            call MatAssemblyBegin(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
            call KSPSetOperators(ksp_solvers(i), A_MAT(i), A_MAT(i), ierr)
            call KSPSetFromOptions(ksp_solvers(i), ierr)
        end if
    end do

    write(*,*) merge(">>> Starting Power Iteration ...", &
                     ">>> Solving for Neutron Flux ...", is_eigenvalue_problem)

    k_eff = 1.0_dp    
    do outer_iter = 1, max_outer_iter
        k_eff_prime = k_eff

        do inner_iter = 1, max_inner_iter
            do i = 1, n_groups
                if (solver_choice == SOLVER_PETSC) then
                    call assemble_multigroup_source_petsc(B_VEC(i), mesh, FE, Quad, materials, X_VEC, X_PRIME_VEC, n_groups, k_eff_prime, i, is_eigenvalue_problem, is_adjoint, S_ext)
                    call KSPSolve(ksp_solvers(i), B_VEC(i), X_VEC(i), ierr)
                else
                    call assemble_multigroup_source_pcg(B_PCG(i), mesh, FE, Quad, materials, X_PCG, X_PRIME_PCG, n_groups, k_eff_prime, i, is_eigenvalue_problem, is_adjoint, S_ext)
                    call CGSolve(A_PCG(i), X_PCG(i), B_PCG(i), preconditioner_choice, max_CG_iter)
                end if
            end do
        end do

        if (.not. is_eigenvalue_problem) exit

        if (solver_choice == SOLVER_PETSC) then
            call calculate_total_production_petsc(k_eff, X_VEC, mesh, FE, Quad, materials, n_groups)
        else
            call calculate_total_production_pcg(k_eff, X_PCG, mesh, FE, Quad, materials, n_groups)
        end if

        do i = 1, n_groups
            if (solver_choice == SOLVER_PETSC) then
                call VecScale(X_VEC(i), 1.0_dp / k_eff, ierr)
                call VecCopy(X_VEC(i), X_PRIME_VEC(i), ierr)
            else
                X_PCG(i)%vec = X_PCG(i)%vec / k_eff
                X_PRIME_PCG(i)%vec = X_PCG(i)%vec
            end if
        end do

        write(*,'(A,I3,A,F12.8)') " Outer Iteration: ", outer_iter, "  k-eff: ", k_eff

        if (abs(k_eff - k_eff_prime) < 1e-7) exit

    end do

    if (solver_choice == SOLVER_PETSC) then
        call export_vtk_petsc("solution_output", FE, mesh, X_VEC, n_groups, 2, .false.)
        call PetscFinalize(ierr)
    else
        call export_vtk_pcg("solution_output", FE, mesh, X_PCG, n_groups, 2, .true.)
    end if

    call system("gprof exe > analysis.txt && python3 -m gprof2dot -f prof analysis.txt | dot -Tsvg -o ../output/profile_map.svg") 

end program fem2d_main