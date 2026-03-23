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
  use m_integrals,  only: integrate_scalar_over_mesh, l2_error_sq_integrand
  use m_petsc
  implicit none

  type(t_material), allocatable :: mats(:)
  type(t_fem_1d)    :: FE
  type(t_quadrature):: quad
  type(t_CSRMatrix), allocatable :: MAT(:)

  integer :: n_groups, gg, nodes_per_elem_plot, n_elem, order, n_quad
  integer :: nnz, iter
  integer, allocatable :: row_ptr(:), col_ind(:), geom, elem_mat_map(:)
  integer, allocatable :: csr_pos(:,:,:)
  real(dp), allocatable :: Nq(:,:), dNq(:,:)
  real(dp), allocatable :: flux_prev(:,:), flux_new(:,:)
  real(dp), allocatable :: K_eff(:)
  real(dp) :: norm_factor, volume_factor
  real(dp) :: keff_den, keff_num, prod_total
  logical :: adjoint
  
  integer :: bc_left, bc_right, ii
  real(dp) :: par_left, par_right
  real(dp) :: x0, x1, t1, t2
  real(dp), allocatable :: VAL_base(:,:)

  call execute_command_line("pkill gnuplot")

  call cpu_time(t1)

  x0     = 0.0_dp
  x1     = 1.0_dp
  adjoint =.false.
  n_groups      = 7

  n_elem = 10
  geom   = GEOM_SLAB
  order  = 1

  bc_left  = BC_NEUMANN
  bc_right = BC_VACUUM
  par_left = 0.0_dp
  par_right= 0.0_dp

  allocate(elem_mat_map(n_elem)); elem_mat_map = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

  nodes_per_elem_plot = 2

  call InitialiseMaterials(mats, n_groups, elem_mat_map, '../MATS.txt', .true.)
  call init_reference_element(FE, order)
  call generate_mesh(FE, x0, x1, n_elem)

  call build_CSR_sparsity_pattern(FE%elements, FE%n_nodes, FE%n_elem, FE%order, nnz, row_ptr, col_ind)

  allocate(MAT(n_groups))
  do gg = 1, n_groups
    allocate(MAT(gg)%row_ptr(FE%n_nodes+1), MAT(gg)%col_ind(nnz), MAT(gg)%val(nnz), MAT(gg)%F(FE%n_nodes), MAT(gg)%u(FE%n_nodes))
    MAT(gg)%row_ptr = row_ptr
    MAT(gg)%col_ind = col_ind
    MAT(gg)%val     = 0.0_dp
    MAT(gg)%F       = 0.0_dp
    MAT(gg)%u       = 1.0_dp
  end do

  allocate(K_eff(10000)); K_eff = 1.0_dp
  allocate(csr_pos(FE%n_elem, FE%nloc, FE%nloc))
  call build_element_csr_positions(FE%n_elem, FE%nloc, FE%elements, row_ptr, col_ind, csr_pos)
  
  n_quad = FE%order + 1
  allocate(quad%Xi(n_quad), quad%W(n_quad))
  call Get1DLineQuad(FE%order+1, quad)
  call precompute_shapes(FE, quad, Nq, dNq)

  allocate(VAL_base(n_groups, nnz))
  do gg = 1, n_groups
    call build_group_matrix(mats, elem_mat_map, FE, quad, Nq, dNq, gg, csr_pos, VAL_base(gg,:), geom)
  end do

  allocate(flux_prev(n_groups, FE%n_nodes), flux_new(n_groups, FE%n_nodes))
  do ii =1,n_elem
  end do

  do iter = 2, 1000

    do gg = 1, n_groups
      flux_prev(gg,:) = MAT(gg)%u(:)
    end do

    call integrate_scalar_over_mesh(mats, elem_mat_map, FE, quad, Nq, dNq, flux_prev, l2_error_sq_integrand, keff_den, adjoint, geom)

    do gg = 1, n_groups
      MAT(gg)%val = VAL_base(gg,:)  ! reuse matrix
      call assemble_group_rhs(mats, elem_mat_map, FE, quad, Nq, dNq, gg, flux_prev, K_eff(iter-1), adjoint, MAT(gg)%F, geom)
      call apply_bcs_csr(MAT(gg)%val, MAT(gg)%col_ind, MAT(gg)%row_ptr, MAT(gg)%F, FE%n_nodes, bc_left, par_left, bc_right, par_right, geom, x1)
      call run_petsc(FE%n_nodes, FE%n_nodes, MAT(gg)%val, MAT(gg)%col_ind, MAT(gg)%row_ptr, MAT(gg)%u, MAT(gg)%F)
    end do

    do gg = 1, n_groups
      flux_new(gg,:) = MAT(gg)%u(:)
    end do

    call integrate_scalar_over_mesh(mats, elem_mat_map, FE, quad, Nq, dNq, flux_new, l2_error_sq_integrand, keff_num, adjoint, geom)
    K_eff(iter) = K_eff(iter-1) * (keff_num / keff_den)

    call integrate_scalar_over_mesh(mats, elem_mat_map, FE, quad, Nq, dNq, flux_new, l2_error_sq_integrand, prod_total, adjoint, geom)

    if (prod_total > 0.0_dp) then
      volume_factor = 1.0_dp
      select case(geom)
          case(GEOM_SPH);  volume_factor = 4.0_dp * PI
          case(GEOM_CYL);  volume_factor = 2.0_dp * PI ! Assumes height H=1cm
          case(GEOM_SLAB); volume_factor = 1.0_dp      ! Assumes area A=1cm^2
      end select

      norm_factor = (1.0_dp / volume_factor) / prod_total
      do gg = 1, n_groups
        MAT(gg)%u = MAT(gg)%u * norm_factor
      end do
    end if

    if (abs(K_eff(iter) - K_eff(iter-1)) < 1.0e-8_dp) exit
  end do

  print*, "K_eff = ", K_eff(iter-1), " Achieved on iteration ", iter

  call cpu_time(t2)
  print*,"Operation Time: ",(t2-t1)

  call write_phi_continuous('phi_continuous.txt', FE, MAT, nodes_per_elem_plot)

call execute_command_line( &
    "gnuplot -e ""set terminal qt size 1100,500 font 'Helvetica,14'; " // &
    "set border lw 1.5 lc rgb '#444444'; set grid lw 1 lc rgb '#dcdcdc'; " // &
    "set xlabel 'x (cm)'; set ylabel 'phi_g'; set key outside; " // &
    "plot 'phi_continuous.txt' index 2 u 1:2 w lp lw 2 lc rgb '#E41A1C' title 'g=1', " // &
    "     'phi_continuous.txt' index 5 u 1:2 w lp lw 2 lc rgb '#377EB8' title 'g=2', " // &
    "     'phi_continuous.txt' index 8 u 1:2 w lp lw 2 lc rgb '#4DAF4A' title 'g=3', " // &
    "     'phi_continuous.txt' index 11 u 1:2 w lp lw 2 lc rgb '#984EA3' title 'g=4', " // &
    "     'phi_continuous.txt' index 14 u 1:2 w lp lw 2 lc rgb '#FF7F00' title 'g=5', " // &
    "     'phi_continuous.txt' index 17 u 1:2 w lp lw 2 lc rgb '#FFFF33' title 'g=6', " // &
    "     'phi_continuous.txt' index 20 u 1:2 w lp lw 2 lc rgb '#A65628' title 'g=7'; " // &
    "pause mouse close""" )
end program fem_1d