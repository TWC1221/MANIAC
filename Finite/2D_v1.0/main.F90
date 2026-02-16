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
    use m_outpVTK
    use m_PCG
    use m_petsc
    
    implicit none
    Mat, allocatable        :: A_MAT(:)    
    Vec, allocatable        :: B_VEC(:)    
    Vec, allocatable        :: X_VEC(:)    
    KSP, allocatable        :: ksp_solvers(:)
    PetscErrorCode          :: ierr

    type(t_mat), allocatable :: A_PCG(:)
    type(t_vec), allocatable :: B_PCG(:)
    type(t_vec), allocatable :: X_PCG(:)

    type(t_mesh)            :: mesh
    type(t_quadrature)      :: Quad, QuadBound
    type(t_finite)          :: FE
    type(t_material), allocatable :: materials(:)
    
    integer                 :: n_groups, i, k, plot_n1d
    integer                 :: n_found
    integer, allocatable    :: vacuum_indices(:), unique_ids(:)
    logical, allocatable    :: is_vacuum_node(:)
    real(dp), allocatable   :: S_ext(:,:)

    integer :: solver_choice, PCG_mode
    integer, parameter :: SOLVER_PCG = 1, SOLVER_PETSC = 2
    integer :: max_iter

    PCG_mode = PCG_PRECON_ILU
    solver_choice = SOLVER_PETSC
    max_iter = 1000

    n_groups = 1
  
    if (solver_choice==SOLVER_PETSC) call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call parse("../VTKmesh.vtk", mesh, .true.)

    call InitializeMaterials(materials, n_groups)

    allocate(unique_ids(mesh%n_elems))
    unique_ids = -1
    n_found = 0
    do i = 1, mesh%n_elems
        if (.not. any(unique_ids == mesh%mats(i))) then
            n_found = n_found + 1
            unique_ids(n_found) = mesh%mats(i)
        end if
    end do

    do i = 1, mesh%n_elems
        do k = 1, n_found
            if (mesh%mats(i) == unique_ids(k)) then
                mesh%mats(i) = k
                exit
            end if
        end do
    end do

    if (FE%order == 3) then
        mesh%nloc = 16
    else if (mesh%nloc == 9) then
        FE%order = 2
    else if (mesh%nloc == 4) then
        FE%order = 1
    end if

    call Get1DLineQuad(FE%order + 1, QuadBound)
    call QuadrilateralQuadrature(Quad, FE%order + 1)
    call InitFE(FE, Quad, QuadBound)
    FE%n_basis = mesh%nloc
    FE%order   = nint(sqrt(real(FE%n_basis))) - 1

    write(*,*) ">>> Starting Global Assembly..."

    if (solver_choice==SOLVER_PETSC) call assemble_petsc_matrix(&
        materials, A_MAT, mesh, FE, Quad, n_groups, mesh%mats)
    if (solver_choice==SOLVER_PCG) call assemble_PCG_matrix(&
        materials, A_PCG, mesh, FE, Quad, n_groups, mesh%mats)

    allocate(S_ext(mesh%n_elems, n_groups))
    do i = 1, mesh%n_elems
        if (mesh%mats(i) == 1) then
            S_ext(i, :) = 10.0_dp
        else
            S_ext(i, :) = 0.0_dp
        end if
    end do

    if (solver_choice==SOLVER_PETSC) call assemble_petsc_source_vec(B_VEC, mesh, FE, Quad, n_groups, S_ext)
    if (solver_choice==SOLVER_PCG) call assemble_PCG_source_vec(B_PCG, mesh, FE, Quad, n_groups, S_ext)

    allocate(is_vacuum_node(mesh%n_nodes)); is_vacuum_node = .false.
    do i = 1, mesh%n_edges
        is_vacuum_node(mesh%edges(i,1)) = .true.
        is_vacuum_node(mesh%edges(i,2)) = .true.
        if (mesh%edges(i,3) > 0) then
            is_vacuum_node(mesh%edges(i,3)) = .true.
        end if
    end do

    allocate(vacuum_indices(count(is_vacuum_node)))
    k = 1
    do i = 1, mesh%n_nodes
        if (is_vacuum_node(i)) then
            vacuum_indices(k) = i - 1
            k = k + 1
        end if
    end do

    if (solver_choice == SOLVER_PETSC) then
        do i = 1, n_groups
            call MatZeroRowsColumns(A_MAT(i), size(vacuum_indices), vacuum_indices, &
                                    1.0_dp, PETSC_NULL_VEC, B_VEC(i), ierr)
            do k = 1, size(vacuum_indices)
                call VecSetValue(B_VEC(i), vacuum_indices(k), 0.0_dp, INSERT_VALUES, ierr)
            end do
            call VecAssemblyBegin(B_VEC(i), ierr)
            call VecAssemblyEnd(B_VEC(i), ierr)
        end do
    else
        do i = 1, n_groups
            call PCG_Apply_Dirichlet(A_PCG(i), B_PCG(i), vacuum_indices)
        end do
    end if

    write(*,*) ">>> Solving for Neutron Flux..."
    if (solver_choice == SOLVER_PETSC) then
        allocate(X_VEC(n_groups))
        allocate(ksp_solvers(n_groups))
    else if (solver_choice == SOLVER_PCG) then
        allocate(X_PCG(n_groups))
        do i = 1, n_groups
            call PCG_VEC_INIT(mesh%n_nodes, X_PCG(i))
        end do
    end if


    do i = 1, n_groups
        if (solver_choice==SOLVER_PETSC) then 
            call VecDuplicate(B_VEC(i), X_VEC(i), ierr)
            call KSPCreate(PETSC_COMM_WORLD, ksp_solvers(i), ierr)
            call KSPSetOperators(ksp_solvers(i), A_MAT(i), A_MAT(i), ierr)
            
            call KSPSetType(ksp_solvers(i), KSPCG, ierr)
            call KSPSetFromOptions(ksp_solvers(i), ierr)
            
            call KSPSolve(ksp_solvers(i), B_VEC(i), X_VEC(i), ierr)
    end if
        else
            call CGSolve(A_PCG(i), X_PCG(i), B_PCG(i), PCG_mode, max_iter)
        end if
        write(*,*) "Group ", i, " converged."
    end do 

    if (solver_choice==SOLVER_PETSC) call export_vtk("solution_output", FE, mesh, X_VEC, n_groups, 20, .true.)
    if (solver_choice==SOLVER_PCG)   call export_vtk_pcg("solution_output", FE, mesh, X_PCG, n_groups, 20, .true.)

    if (solver_choice==SOLVER_PETSC) then
    do i = 1, n_groups
        call MatDestroy(A_MAT(i), ierr)
        call VecDestroy(B_VEC(i), ierr)
        call VecDestroy(X_VEC(i), ierr)
        call KSPDestroy(ksp_solvers(i), ierr)
    end do
    call PetscFinalize(ierr)
    end if

    write(*,*) ">>> Simulation Finished Successfully."

end program fem2d_main