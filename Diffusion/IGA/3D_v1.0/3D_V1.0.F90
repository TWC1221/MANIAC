#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

program fem2d_main
    use m_constants
    use m_types
    use m_asmg
    use m_quadrature
    use m_basis
    use m_material
    use m_boundaries
    use m_outpVTK
    use m_PCG
    use m_petsc

    implicit none

    Mat, allocatable        :: A_MAT(:)    
    Mat, allocatable        :: MAT_F_PETSC(:,:), MAT_S_PETSC(:,:)
    Vec, allocatable        :: B_VEC(:), X_VEC(:), X_PRIME_VEC(:)    
    Vec, allocatable        :: FixedSrc_PETSC(:), PROD_VEC(:)
    KSP, allocatable        :: ksp_solvers(:)
    PetscErrorCode          :: ierr

    type(t_PCG_CSR), allocatable :: A_PCG(:)
    type(t_PCG_CSR), allocatable :: MAT_F_PCG(:,:), MAT_S_PCG(:,:)
    type(t_vec), allocatable :: B_PCG(:), X_PCG(:), X_PRIME_PCG(:), FixedSrc_PCG(:)

    type(t_mesh)            :: mesh
    type(t_bc_config)       :: bc_config(2)
    type(t_quadrature)      :: Quad, QuadBound
    type(t_finite)          :: FE
    type(t_material), allocatable :: materials(:)
    
    integer                 :: n_groups, i, k, outer_iter
    real(dp)                :: max_phi_change, group_diff
    integer                 :: solver_choice, preconditioner_choice
    integer                 :: max_outer_iter, max_CG_iter
    logical                 :: is_eigenvalue_problem, is_adjoint
    real(dp)                :: k_eff, k_eff_prime, total_production, mesh_volume
    real(dp), allocatable   :: mat_volumes(:)
    character(len=32)       :: InputMesh

    InputMesh = "../input/rod_test.asmg"
    is_eigenvalue_problem   = .true.
    is_adjoint              = .false.
    n_groups               = 7

    solver_choice          = SOLVER_KSP_CG
    preconditioner_choice  = PRECON_GAMG
    max_outer_iter         = 500          
    max_CG_iter            = 100000          

    call read_asmg_mesh(InputMesh, mesh)
    FE%n_basis = mesh%nloc
    FE%order   = mesh%order

    call Get1DLineQuad(QuadBound, FE%order + 1)
    call HexahedralQuadrature(Quad, QuadBound, FE%order + 1)
    call QuadrilateralQuadrature(QuadBound, FE%order + 1)

    call InitialiseBasis(FE)
    call InitialiseMaterials(materials, mesh, n_groups, "../input/MATS.txt", printout = .false.)

    call calculate_mesh_volume(mesh_volume, mat_volumes, mesh, FE, Quad)
    write(*,'(A, F15.6)') " >>> Total Mesh Volume:    ", mesh_volume
    do i = 1, size(mat_volumes)
        if (mat_volumes(i) > 0.0_dp) then
            write(*,'(A, I3, A, F15.6)') " >>> Volume Material ", i, ": ", mat_volumes(i)
        end if
    end do

    call InitialiseBoundaries(bc_config(1), 3, BC_VACUUM, 0.0_dp)
    call InitialiseBoundaries(bc_config(2), 2, BC_REFLECTIVE, 0.0_dp)

    if (solver_choice /= SOLVER_PCG) then
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        allocate(A_MAT(n_groups), B_VEC(n_groups), X_VEC(n_groups), X_PRIME_VEC(n_groups), PROD_VEC(n_groups), FixedSrc_PETSC(n_groups), ksp_solvers(n_groups))
        do i = 1, n_groups
            call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, B_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), X_VEC(i), ierr)
            call VecDuplicate(B_VEC(i), X_PRIME_VEC(i), ierr)
            call VecSet(X_VEC(i), merge(1.0_dp, 0.0_dp, is_eigenvalue_problem), ierr)
            call VecSet(X_PRIME_VEC(i), merge(1.0_dp, 0.0_dp, is_eigenvalue_problem), ierr)  
            call KSPCreate(PETSC_COMM_SELF, ksp_solvers(i), ierr)
        end do
        call ASSEMBLEMultigroupMAT_PETSc(A_MAT, MAT_F_PETSC, MAT_S_PETSC, PROD_VEC, FixedSrc_PETSC, mesh, FE, Quad, materials, n_groups, is_adjoint)
    else 
        allocate(A_PCG(n_groups), B_PCG(n_groups), X_PCG(n_groups), X_PRIME_PCG(n_groups))
        do i = 1, n_groups
            call PCG_VEC_INIT(mesh%n_nodes, X_PCG(i))
            call PCG_VEC_INIT(mesh%n_nodes, B_PCG(i))
            call PCG_VEC_INIT(mesh%n_nodes, X_PRIME_PCG(i))
            X_PCG(i)%vec = merge(1.0_dp, 0.0_dp, is_eigenvalue_problem)
            X_PRIME_PCG(i)%vec = merge(1.0_dp, 0.0_dp, is_eigenvalue_problem)
        end do
         write(*,*) ">>> Assembling Global Matrix ..."
        call assemble_PCG_matrix(materials, A_PCG, mesh, FE, Quad, n_groups, mesh%mats)
        call assemble_source_matrices_pcg(MAT_F_PCG, MAT_S_PCG, FixedSrc_PCG, mesh, FE, Quad, materials, n_groups, is_adjoint)
    end if

    ! Apply Boundary Conditions
    do i = 1, n_groups
        do k = 1, size(bc_config)
            if (solver_choice /= SOLVER_PCG) then
                call apply_bcs(mesh, FE, QuadBound, bc_config(k), A_petsc=A_MAT(i))
            else
                call apply_bcs(mesh, FE, Quadbound, bc_config(k), A_pcg=A_PCG(i))
            end if
        end do

            ! Finalize PETSc Assembly after domain and BCs are applied
            if (solver_choice /= SOLVER_PCG) then
                call MatAssemblyBegin(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyEnd(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
                do k = 1, n_groups
                    call MatAssemblyBegin(MAT_F_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
                    call MatAssemblyEnd(MAT_F_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
                    call MatAssemblyBegin(MAT_S_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
                    call MatAssemblyEnd(MAT_S_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
                end do
            end if

        if (solver_choice /= SOLVER_PCG) then
            call SetupKSP_PETSc(ksp_solvers(i), A_MAT(i), solver_choice, preconditioner_choice)
        end if
    end do

    write(*,*) merge(">>> Starting Power Iteration ...", &
                     ">>> Solving for Neutron Flux ...", is_eigenvalue_problem)

    k_eff = 1.0_dp     
    do outer_iter = 1, max_outer_iter
        k_eff_prime = k_eff
        max_phi_change = 0.0_dp

        do i = 1, n_groups
            if (solver_choice /= SOLVER_PCG) then
                call EVALMultigroup_PETSc(B_VEC(i), MAT_F_PETSC, MAT_S_PETSC, FixedSrc_PETSC, &
                                                      X_VEC, X_PRIME_VEC, n_groups, k_eff_prime, i, is_eigenvalue_problem)
            else
                call assemble_multigroup_source_pcg(B_PCG(i), MAT_F_PCG, MAT_S_PCG, FixedSrc_PCG, &
                                                    X_PCG, X_PRIME_PCG, n_groups, k_eff_prime, i, is_eigenvalue_problem)
            end if
        end do

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

            write(*,'(A,I3,A,F12.8,A,F12.8)') " Power Iteration: ", outer_iter, " k-eff: ", k_eff, "        Delta Phi: ", max_phi_change
            if (abs(k_eff - k_eff_prime) < 1e-5) exit
        else
            write(*,'(A,I3,A,ES12.4)') " Source Iteration: ", outer_iter, " Max Change: ", max_phi_change
            if (max_phi_change < 1e-8) exit
        end if 
    end do

    if (solver_choice /= SOLVER_PCG) then
        do i = 1, n_groups
            call MatAssemblyBegin(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(A_MAT(i), MAT_FINAL_ASSEMBLY, ierr)
            do k = 1, n_groups
                call MatAssemblyBegin(MAT_F_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyEnd(MAT_F_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyBegin(MAT_S_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyEnd(MAT_S_PETSC(i, k), MAT_FINAL_ASSEMBLY, ierr)
            end do
        end do
    end if

    if (solver_choice /= SOLVER_PCG) then
        call export_vtk_petsc("../output/"//derive_case_nametag(InputMesh), FE, mesh, X_VEC, n_groups, 10, .false.)
        call PetscFinalize(ierr)
    else
        call export_vtk_pcg("../output/"//derive_case_nametag(InputMesh), FE, mesh, X_PCG, n_groups, 10, .false.)
    end if
    
    call system("gprof exe > analysis.txt && python3 -m gprof2dot -f prof analysis.txt | dot -Tsvg -o ../output/profile_map"//derive_case_nametag(trim(InputMesh))//".svg") 

end program fem2d_main
