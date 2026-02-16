program fem_1d
  use m_constants
  use m_materials
  use m_fem
  use m_quadrature
  use m_CSR_types
  use m_csr_utils
  use m_assembly
  use m_output
  use m_boundaries
  use m_PCG_solver
  use m_integrals,  only: integrate_scalar_over_mesh, fission_integrand, production_integrand
  use m_petsc
  implicit none

  type(t_material)  :: mats
  type(t_fem_1d)    :: fem
  type(t_quadrature):: quad
  type(t_CSRMatrix), allocatable :: MAT(:)

  integer :: G, gg, nodes_per_elem_plot, n_elem, order, n_quad
  integer :: nnz, iter
  integer, allocatable :: row_ptr(:), col_ind(:)
  integer, allocatable :: csr_pos(:,:,:)
  real(dp), allocatable :: Nq(:,:), dNq(:,:)
  real(dp), allocatable :: flux_prev(:,:), flux_new(:,:)
  real(dp), allocatable :: K_eff(:)
  real(dp) :: keff_den, keff_num, prod_total
  logical :: adjoint
  integer, parameter :: GEOM_SLAB=0, GEOM_CYL=1, GEOM_SPH=2
  integer, parameter :: BC_NONE=0, BC_DIRICHLET=1, BC_NEUMANN=2, BC_VACUUM=3, BC_ALBEDO=4, BC_PERIODIC=5
  integer :: bc_left, bc_right, ii
  real(dp) :: par_left, par_right
  real(dp) :: x0, x1, t1, t2
  real(dp), allocatable :: VAL_base(:,:)

  call execute_command_line("pkill gnuplot")

  call cpu_time(t1)

  ! Problem setup
  adjoint =.false.
  G      = 3
  call init_material_const(mats, G, GEOM_SPH)

  bc_left  = 2  ! BC_NEUMANN
  bc_right = 3  ! BC_VACUUM
  par_left = 0.0_dp
  par_right= 0.0_dp

  x0     = 0.0_dp
  x1     = 10.0_dp
  order = 2
  n_elem = 4

  nodes_per_elem_plot = 6

  call init_reference_element(fem, order)
  call generate_mesh(fem, x0, x1, n_elem)
  call build_CSR_sparsity_pattern(fem%elements, fem%n_nodes, fem%n_elem, fem%order, nnz, row_ptr, col_ind)

  allocate(MAT(G))
  do gg = 1, G
    allocate(MAT(gg)%row_ptr(fem%n_nodes+1), MAT(gg)%col_ind(nnz), MAT(gg)%val(nnz), MAT(gg)%F(fem%n_nodes), MAT(gg)%u(fem%n_nodes))
    MAT(gg)%row_ptr = row_ptr
    MAT(gg)%col_ind = col_ind
    MAT(gg)%val     = 0.0_dp
    MAT(gg)%F       = 0.0_dp
    MAT(gg)%u       = 1.0_dp
  end do

  allocate(K_eff(10000)); K_eff = 1.0_dp
  allocate(csr_pos(fem%n_elem, fem%nloc, fem%nloc))
  call build_element_csr_positions(fem%n_elem, fem%nloc, fem%elements, row_ptr, col_ind, csr_pos)
  
  n_quad = fem%order + 1
  allocate(quad%Xi(n_quad), quad%W(n_quad))
  call Get1DLineQuad(fem%order+1, quad)
  call precompute_shapes(fem, quad, Nq, dNq)

  allocate(VAL_base(G, nnz))
  do gg = 1, G
    call build_group_matrix_once(mats, fem, quad, Nq, dNq, gg, csr_pos, VAL_base(gg,:))
  end do

  allocate(flux_prev(G, fem%n_nodes), flux_new(G, fem%n_nodes))
  do ii =1,n_elem
  print*,fem%elements(ii,:)
  end do

  do iter = 2, 1000

    do gg = 1, G
      flux_prev(gg,:) = MAT(gg)%u(:)
    end do

    call integrate_scalar_over_mesh(mats, fem, quad, Nq, dNq, flux_prev, fission_integrand, keff_den, adjoint)

    do gg = 1, G
      MAT(gg)%val = VAL_base(gg,:)  ! reuse matrix
      call assemble_group_rhs_only(mats, fem, quad, Nq, dNq, gg, flux_prev, K_eff(iter-1), adjoint, MAT(gg)%F)
      call apply_bcs_csr(MAT(gg)%val, MAT(gg)%col_ind, MAT(gg)%row_ptr, MAT(gg)%F, fem%n_nodes, bc_left, par_left, bc_right, par_right, mats%geom_m, x1)
      call PCG_algorithm(MAT(gg)%val, MAT(gg)%col_ind, MAT(gg)%row_ptr, MAT(gg)%u, MAT(gg)%F, 0, 20000)
    end do

    do gg = 1, G
      flux_new(gg,:) = MAT(gg)%u(:)
    end do

    call integrate_scalar_over_mesh(mats, fem, quad, Nq, dNq, flux_new, fission_integrand, keff_num, adjoint)
    K_eff(iter) = K_eff(iter-1) * (keff_num / keff_den)

    call integrate_scalar_over_mesh(mats, fem, quad, Nq, dNq, flux_new, production_integrand, prod_total, adjoint)

    if (prod_total > 0.0_dp) then
      do gg = 1, G
        MAT(gg)%u = MAT(gg)%u * K_eff(iter) / prod_total
      end do
    end if

    if (abs(K_eff(iter) - K_eff(iter-1)) < 1.0e-8_dp) exit
  end do

  print*, K_eff(iter-1), iter

  call cpu_time(t2)
  print*,(t2-t1)

  call write_phi_continuous('phi_continuous.txt', fem, MAT, nodes_per_elem_plot)

  call execute_command_line( &
    "gnuplot -e ""set terminal qt size 1100,500 font 'Helvetica,14'; " // &
    "set border lw 1.5 lc rgb '#444444'; set grid lw 1 lc rgb '#dcdcdc'; " // &
    "set xlabel 'x (cm)'; set ylabel 'ϕ_{g}'; set key outside; " // &
    "plot 'phi_continuous.txt' index 2 u 1:2 w lp lw 2 lc rgb '#E41A1C' title 'g=1 interp', " // &
    "     'phi_continuous.txt' index 5 u 1:2 w lp lw 2 lc rgb '#377EB8' title 'g=2 interp', " // &
    "     'phi_continuous.txt' index 8 u 1:2 w lp lw 2 lc rgb '#4DAF4A' title 'g=3 interp'; " // &
    "pause mouse close""" )
end program fem_1d