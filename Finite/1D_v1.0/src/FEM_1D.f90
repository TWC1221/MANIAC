program fem_1d
  use m_constants
  use m_finite_element
  use m_quadrature
  use m_CSR_types
  use m_boundaries
  use m_PCG_solver
  use m_petsc
  implicit none

  type(t_quadrature) :: quad
  type(t_CSRMatrix), allocatable :: MAT(:)

  ! FEM objects
  type(t_fem_1d) :: fe
  real(dp), allocatable :: x_nodes(:)
  integer,  allocatable :: elements(:,:)

  ! Local variables
  integer :: G, gg, ee, ii, jj, kk, q
  real(dp) :: xi, wq, xq, J
  real(dp) :: Dq, Sr
  integer :: jg
  real(dp) :: keff_den, keff_num
  real(dp) :: scatter_in_q, fission_prev, fission_new, prod_total
  real(dp) :: phi_prev_jg, phi_new_jg
  integer  :: n_elem, order, n_nodes, nodes_per_elem_plot
  real(dp) :: x0, x1

  ! CSR VARS
  integer :: idx, jdx, nnz
  integer, allocatable :: row_ptr(:), col_ind(:)
  real(dp), allocatable :: val(:), Kel(:,:), Fel(:), xe(:)

  integer :: iter, nloc
  real(dp), allocatable :: K_eff(:), u_prev(:,:)
  real(dp), allocatable :: Ni(:), dNi(:), ue(:)

  ! BC
  integer, parameter :: BC_NONE=0, BC_DIRICHLET=1, BC_NEUMANN=2, BC_VACUUM=3, BC_ALBEDO=4, BC_PERIODIC=5
  integer :: bc_left, bc_right
  real(dp) :: par_left, par_right, t1, t2

  ! Geometry
  integer, parameter :: GEOM_SLAB=0, GEOM_CYL=1, GEOM_SPH=2
  integer :: geom_m
  real(dp) :: geo_w
  real(dp), parameter :: TINY_R = 1.0e-14_dp

  logical :: adjoint

  call execute_command_line("pkill gnuplot")

  call cpu_time(t1)

  adjoint = .false.
  geom_m = GEOM_SPH
  bc_left  = BC_NEUMANN;  par_left  = 0.0_dp
  bc_right = BC_VACUUM;  par_right = 0.0_dp

  x0     = 0.0_dp
  x1     = 10.0_dp
  n_elem = 4
  order  = 2
  G      = 3

  nodes_per_elem_plot = 6

  call init_reference_element_1d(fe, order)
  call generate_mesh_1d(fe, x0, x1, n_elem, x_nodes, elements)
  n_nodes = size(x_nodes)
  call build_CSR_sparsity_pattern(elements, n_nodes, n_elem, order, nnz, row_ptr, col_ind)

  allocate(MAT(G), val(nnz))
  allocate(Quad%Xi(order+1))
  allocate(Quad%W(order+1))
  allocate(Kel(order+1,order+1), Fel(order+1), xe(order+1), K_eff(10000))

  do gg = 1, G
    allocate(MAT(gg)%row_ptr(n_nodes+1), MAT(gg)%col_ind(nnz), MAT(gg)%val(nnz), MAT(gg)%F(n_nodes), MAT(gg)%u(n_nodes))
    MAT(gg)%row_ptr = row_ptr
    MAT(gg)%col_ind = col_ind
    MAT(gg)%val     = 0.0_dp
    MAT(gg)%F       = 0.0_dp
    MAT(gg)%u       = 1.0_dp
  end do
  K_eff = 1.0_dp

  nloc = order + 1
  allocate(u_prev(G, n_nodes), Ni(nloc), dNi(nloc), ue(nloc))
  do gg = 1, G
    u_prev(gg,:) = MAT(gg)%u(:)
  end do 

  do iter = 2, 1000

    call Get1DLineQuad(order+1, quad)

    keff_den = 0.0_dp
    do ee = 1, n_elem
      xe(:) = x_nodes(elements(ee,1:nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)

        do ii = 1, nloc
          Ni(ii)  = lagrange_shape_1d(fe, ii, xi)
          dNi(ii) = lagrange_shape_derivative_1d(fe, ii, xi)
        end do
        
        xq = dot_product(Ni,  xe)
        J  = dot_product(dNi, xe)
        geo_w = max(xq, TINY_R)**geom_m

        fission_prev = 0.0_dp
        do jg = 1, G
          ue(:) = u_prev(jg, elements(ee,1:nloc))
          phi_prev_jg = dot_product(ue, Ni)
          if (.not. adjoint) fission_prev = fission_prev + Sigmanuf_of_x(xq, jg) * phi_prev_jg
          if (adjoint) fission_prev = fission_prev + chi_of_x(xq, jg) * phi_prev_jg
        end do

        keff_den = keff_den + wq * geo_w * fission_prev * J
      end do
    end do

    ! --- Group solves ---
    do gg = 1, G
      MAT(gg)%val = 0.0_dp
      MAT(gg)%F   = 0.0_dp

      do ee = 1, n_elem
        Kel = 0.0_dp
        Fel = 0.0_dp
        xe(:) = x_nodes(elements(ee,1:nloc))

        do q = 1, quad%NoPoints
          xi = quad%Xi(q); wq = quad%W(q)

          do ii = 1, nloc
            Ni(ii)  = lagrange_shape_1d(fe, ii, xi)
            dNi(ii) = lagrange_shape_derivative_1d(fe, ii, xi)
          end do

          xq = dot_product(Ni,  xe)
          J  = dot_product(dNi, xe)
          geo_w = max(xq, TINY_R)**geom_m

          Dq = D_of_x(xq, gg, G)
          Sr = SigmaR_of_x(xq, gg, G)

          scatter_in_q = 0.0_dp
          fission_prev = 0.0_dp
          do jg = 1, G

            ue(:) = u_prev(jg, elements(ee,1:nloc))
            phi_prev_jg = dot_product(ue, Ni)

            if (jg /= gg) then
              if (.not. adjoint) then
                scatter_in_q = scatter_in_q + SigmaS_of_x(xq, jg, gg) * phi_prev_jg
              else
                scatter_in_q = scatter_in_q + SigmaS_of_x(xq, gg, jg) * phi_prev_jg
              end if
            end if

            if (.not. adjoint) fission_prev = fission_prev + Sigmanuf_of_x(xq, jg) * phi_prev_jg
            if (adjoint) fission_prev = fission_prev + chi_of_x(xq, jg) * phi_prev_jg

          end do

          do ii = 1, nloc
            do jj = 1, nloc
              Kel(ii,jj) = Kel(ii,jj) + wq * geo_w * ( Dq * (dNi(ii)*dNi(jj)) / J + Sr * Ni(ii) * Ni(jj) * J )
            end do
            if (.not. adjoint) Fel(ii) = Fel(ii) + wq * geo_w * ( scatter_in_q + chi_of_x(xq, gg) * fission_prev / K_eff(iter-1) ) * Ni(ii) * J
            if (adjoint) Fel(ii) = Fel(ii) + wq * geo_w * ( scatter_in_q + Sigmanuf_of_x(xq, gg) * fission_prev / K_eff(iter-1) ) * Ni(ii) * J
          end do
        end do

        do ii = 1, nloc
          MAT(gg)%F(elements(ee,ii)) = MAT(gg)%F(elements(ee,ii)) + Fel(ii)
        end do
        do ii = 1, nloc
          idx = elements(ee,ii)
          do jj = 1, nloc
            jdx = elements(ee,jj)
            do kk = MAT(gg)%row_ptr(idx), MAT(gg)%row_ptr(idx+1)-1
              if (MAT(gg)%col_ind(kk) == jdx) then
                MAT(gg)%val(kk) = MAT(gg)%val(kk) + Kel(ii,jj)
                exit
              end if
            end do
          end do
        end do
      end do  ! elements

      call apply_bcs_csr(MAT(gg)%val, MAT(gg)%col_ind, MAT(gg)%row_ptr, MAT(gg)%F, n_nodes, bc_left, par_left, bc_right, par_right, geom_m, x1)
      call PCG_algorithm(MAT(gg)%val, MAT(gg)%col_ind, MAT(gg)%row_ptr, MAT(gg)%u, MAT(gg)%F, 0, 20000)
      !call run_petsc(n_nodes, n_nodes, MAT(gg)%val, MAT(gg)%col_ind, MAT(gg)%row_ptr, MAT(gg)%u, MAT(gg)%F)
    end do  ! gg

    keff_num = 0.0_dp
    do ee = 1, n_elem
      xe(:) = x_nodes(elements(ee,1:nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)

        do ii = 1, nloc
          Ni(ii)  = lagrange_shape_1d(fe, ii, xi)
          dNi(ii) = lagrange_shape_derivative_1d(fe, ii, xi)
        end do
        xq = dot_product(Ni,  xe)
        J  = dot_product(dNi, xe)
        geo_w = max(xq, TINY_R)**geom_m

        fission_new = 0.0_dp
        do jg = 1, G
          ue(:) = MAT(jg)%u(elements(ee,1:nloc))
          phi_new_jg = dot_product(ue, Ni)
          if (.not. adjoint) fission_new = fission_new + Sigmanuf_of_x(xq, jg) * phi_new_jg
          if (adjoint) fission_new = fission_new + chi_of_x(xq, jg) * phi_new_jg
        end do

        keff_num = keff_num + wq * geo_w * fission_new * J
      end do
    end do

    K_eff(iter) = K_eff(iter-1) * (keff_num / keff_den)

    prod_total = 0.0_dp
    do ee = 1, n_elem
      xe(:) = x_nodes(elements(ee,1:nloc))
      do q = 1, quad%NoPoints
        xi = quad%Xi(q); wq = quad%W(q)

        do ii = 1, nloc
          Ni(ii)  = lagrange_shape_1d(fe, ii, xi)
          dNi(ii) = lagrange_shape_derivative_1d(fe, ii, xi)
        end do
        xq = dot_product(Ni,  xe)
        J  = dot_product(dNi, xe)
        geo_w = max(xq, TINY_R)**geom_m

        do jg = 1, G
          ue(:) = MAT(jg)%u(elements(ee,1:nloc))
          prod_total = prod_total + wq * geo_w * Sigmanuf_of_x(xq, jg) * dot_product(ue, Ni) * J
        end do
      end do
    end do

    if (prod_total > 0.0_dp) then
      do gg = 1, G
        MAT(gg)%u = MAT(gg)%u * K_eff(iter) / prod_total
      end do
    end if

    do gg = 1, G
      u_prev(gg,:) = MAT(gg)%u(:)
    end do

    if (abs(K_eff(iter) - K_eff(iter-1)) < 1.0e-8_dp) exit
  end do  ! iter

  print*, K_eff(iter-1), iter

  call cpu_time(t2)
  print*,(t2-t1)
  call write_phi_continuous('phi_continuous.txt', fe, elements, x_nodes, MAT, nodes_per_elem_plot)

  call execute_command_line( &
    "gnuplot -e ""set terminal qt size 1100,500 font 'Helvetica,14'; " // &
    "set border lw 1.5 lc rgb '#444444'; set grid lw 1 lc rgb '#dcdcdc'; " // &
    "set xlabel 'x (cm)'; set ylabel 'phi_g(x)'; set key outside; " // &
    "plot 'phi_continuous.txt' index 2 u 1:2 w lp lw 2 lc rgb '#E41A1C' title 'g=1 interp', " // &
    "     'phi_continuous.txt' index 5 u 1:2 w lp lw 2 lc rgb '#377EB8' title 'g=2 interp', " // &
    "     'phi_continuous.txt' index 8 u 1:2 w lp lw 2 lc rgb '#4DAF4A' title 'g=3 interp'; " // &
    "pause mouse close""" )

contains

subroutine build_CSR_sparsity_pattern(elements, n_nodes, n_elem, order, nnz, row_ptr, col_ind)
  implicit none

  integer, intent(in)  :: n_nodes, n_elem, order
  integer, intent(in)  :: elements(n_elem, order+1)

  integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)

  integer :: idx, jdx, ee, ii, jj
  integer :: nnz, kk
  integer, allocatable :: marker(:)

  allocate(row_ptr(n_nodes+1), marker(n_nodes))

  kk = 1
  nnz = 0
  marker = 0
  row_ptr(1) = 1

  do idx = 1, n_nodes
    do ee = 1, n_elem
      do ii = 1, order+1
        if (elements(ee,ii) == idx) then
          do jj = 1, order+1
            jdx = elements(ee,jj)
            if (marker(jdx) == 0) then
              marker(jdx) = 1
              nnz = nnz + 1
            end if
          end do
        end if
      end do
    end do
    row_ptr(idx+1) = nnz + 1
    marker = 0
  end do

  allocate(col_ind(nnz))
  marker = 0

  do idx = 1, n_nodes
    do ee = 1, n_elem
      do ii = 1, order+1
        if (elements(ee,ii) == idx) then
          do jj = 1, order+1
            jdx = elements(ee,jj)
            if (marker(jdx) == 0) then
              marker(jdx) = 1
              col_ind(kk) = jdx
              kk = kk + 1
            end if
          end do
        end if
      end do
    end do
    marker = 0
  end do

  deallocate(marker)

end subroutine build_csr_sparsity_pattern

end program FEM_1D

