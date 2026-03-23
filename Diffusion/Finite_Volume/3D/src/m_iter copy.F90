
!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  Conduct 2D NDE eigenvalue problems in xy, rz and rth geometries      -!  
!  Input data defining problem geometry, discretization and meshing     -!
!  Input data defining neutron/material interaction                     -!
!  Call functions to assemble & solve coefficient matrix & output VTK   -!
!    - assemble_XX_diffusion_matrix                                     -!
!    - PCG_algorithm                                                    -!
!    - outpVTK                                                          -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 12/11/25     T. Charlton    Implemented 2D matrix                     -!
! 15/11/25     T. Charlton    Implemented rth matrix                    -!
! 17/11/25     T. Charlton    Implemented rz matrix                     -!
! 27/11/25     T. Charlton    Implemented 3D matrix                     -!
! 07/12/25     T. Charlton    Implemented 3D multigroup                 -!
! 10/12/25     T. Charlton    Implemented 3D current calcs              -! 
!------------------------------------------------------------------------! 

module m_iter
  use CSR_types
  use m_assembly
  use m_PCG_solver
  use m_outpVTK
  use m_wall_current
  use m_kinetic
  use m_subcritical
  use m_current
  use m_petsc
  implicit none
  contains

  subroutine iter_xyz_multi(PCG_mode, max_iter, adjoint, fixed_source, subcritical, kinetics)
    use BC_types
    use mesh_types
    implicit none

    type(CSRMatrix), allocatable :: A_all(:)
    type(MeshGrid) :: mesh
    type(BoundarySet) :: BCs
    type(MatMatrix) :: MATs

    logical, intent(in), optional :: adjoint, fixed_source, subcritical, kinetics
    integer, intent(in) :: PCG_mode, max_iter
    integer :: ii, jj, kk, G, gg

    ! fields and cross-sections
    real(8), allocatable :: phi(:,:), phi_prime(:,:), src_ext(:,:), phi_ptr(:)
    real(8), allocatable :: D(:), Sigma_t(:), Sigma_a(:), Sigma_r(:), nu_Sigma_f(:)

    real(8), allocatable :: Sigma_s(:,:), Sigma_s_L(:,:), Sigma_s_U(:,:)
    real(8), allocatable :: Sigma_s_L_sum(:,:), Sigma_s_U_sum(:,:)
    real(8), allocatable :: S_f(:), S_f_iter(:), K_eff(:), chi(:), fission_rate(:)

    real(8), allocatable :: Jx_min(:,:), Jx_max(:,:)
    real(8), allocatable :: Jy_min(:,:), Jy_max(:,:)
    real(8), allocatable :: Jz_min(:,:), Jz_max(:,:)
    real(8), allocatable :: Jx(:,:), Jy(:,:), Jz(:,:)

    real(8) :: P, Abso, L, norm
    logical :: converged_inner, use_adjoint, use_fixed_source, use_subcritical, use_kinetics

    use_adjoint =.false.
    use_fixed_source=.false.
    use_subcritical=.false.
    use_kinetics =.false.

    if (present(adjoint)) use_adjoint = adjoint
    if (present(fixed_source)) use_fixed_source = fixed_source
    if (present(subcritical)) use_subcritical = subcritical
    if (present(kinetics)) use_kinetics = kinetics

    ! probelm size
    call mesh%init('../mesh_inputs.txt') 
    call read_material_file('../mat_1.1%_U235.txt', MATs)

    G = MATs%G
    Sigma_s = MATs%Sigma_s
    Sigma_a = MATs%Sigma_a
    nu_Sigma_f = MATs%nu_Sigma_f
    chi = MATs%chi

    allocate(A_all(G), phi(G, mesh%N), phi_prime(G, mesh%N), phi_ptr(mesh%N))
    allocate(D(G), Sigma_t(G), Sigma_r(G), src_ext(G,mesh%N))

    allocate(Sigma_s_L(G,G), Sigma_s_U(G,G))
    allocate(Sigma_s_L_sum(G,mesh%N), Sigma_s_U_sum(G,mesh%N))
    allocate(S_f(mesh%N), S_f_iter(mesh%N), K_eff(max_iter+1), fission_rate(mesh%N))

    allocate(Jx(G,mesh%N), Jy(G,mesh%N), Jz(G,mesh%N))
    allocate(Jx_min(G, mesh%ny*mesh%nz), Jx_max(G, mesh%ny*mesh%nz))
    allocate(Jy_min(G, mesh%nx*mesh%nz), Jy_max(G, mesh%nx*mesh%nz))
    allocate(Jz_min(G, mesh%nx*mesh%ny), Jz_max(G, mesh%nx*mesh%ny))

    ! Initialize all faces as Vacuum
    call init_boundary_set(BCs, default_type = 4)  ! 4 = Vacuum

    ! set_bc_dirichlet(BCs, face_ID, value), set_bc_neumann(BCs, face_ID, value), set_bc_vacuum(BCs, face_ID), set_bc_albedo(BCs, face_ID, alpha), periodic at call ASSEMBLY
    call set_bc_dirichlet(BCs, FACE_XMIN, 0.0d0)
    call set_bc_dirichlet(BCs, FACE_XMAX, 0.0d0)
    call set_bc_neumann(BCs, FACE_YMIN, 0.0d0)
    call set_bc_neumann(BCs, FACE_YMAX, 0.0d0)
    call set_bc_neumann(BCs, FACE_ZMIN, 0.0d0)
    call set_bc_neumann(BCs, FACE_ZMAX, 0.0d0)

    ! initial flux guess
    phi = 1.0d-6
    phi_prime = 1.0d-6
    S_f = 1.0d0
    K_eff(1) = 1.0d0

    do gg = 1, G
      Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s(gg,1:G))
      Sigma_r(gg) = Sigma_t(gg) - Sigma_s(gg,gg)
      D(gg) = 1.0d0 / (3.0d0 * Sigma_t(gg))
      call assemble_3D_diffusion_matrix_g(D(gg), Sigma_r(gg), A_all(gg), mesh, BCs, use_adjoint, periodic_x = .false., periodic_y = .false., periodic_z = .false. )
    end do

    if (use_adjoint) Sigma_s = transpose(Sigma_s)

    ! split into L (downscatter to lower groups) and U (upscatter to higher groups)
    Sigma_s_L = 0.0d0
    Sigma_s_U = 0.0d0

    do ii = 1, G
      do jj = 1, G
        if (ii > jj) Sigma_s_L(ii,jj) = Sigma_s(ii,jj)  ! i>j: downscatter (to lower group index)
        if (ii < jj) Sigma_s_U(ii,jj) = Sigma_s(ii,jj)  ! i<j: upscatter
      end do
    end do

    if (use_adjoint .and. use_subcritical) then
      write(*,'(A)') char(27)//"[1;34m"//"===== NEUTRON IMPORTANCE SOLVE ====="//char(27)//"[0m"
      do ii = 1, max_iter
        phi_prime = phi
        do gg = 1, G

          Sigma_s_L_sum = 0.0d0
          Sigma_s_U_sum = 0.0d0

          do jj = 1, G
            Sigma_s_L_sum(jj,:) = Sigma_s_L(jj,gg) * phi(jj,:)
            Sigma_s_U_sum(jj,:) = Sigma_s_U(jj,gg) * phi_prime(jj,:)
          end do

          S_f = (sum(Sigma_s_L_sum, dim=1) + sum(Sigma_s_U_sum, dim=1) + nu_Sigma_f(gg)) * mesh%dV 

          phi_ptr = phi(gg,:)
          call PCG_algorithm(A_all(gg)%aa, A_all(gg)%ja, A_all(gg)%ia, phi_ptr, S_f, PCG_mode, max_iter)
          phi(gg,:) = phi_ptr

        end do

        if (mod(ii, 10) == 0) then
            print *, ( sum(abs(phi(gg,:) - phi_prime(gg,:))), gg = 1, G )
        end if

        converged_inner =.true.
        do gg = 1, G
          if (sum(abs(phi(gg,:) - phi_prime(gg,:))) > 1.0d-6) then
            converged_inner =.false.
          end if
        end do
        if (converged_inner) exit
        
      end do

      call outpVTK_xyz_multi('G', phi, mesh, G, use_adjoint)

      write(*,'(A)') char(27)//"[1;34m"//"===== NEUTRON IMPORTANCE EVALUATED ====="//char(27)//"[0m"
    end if

    if (use_fixed_source) then
      print *, ""
      write(*,'(A)') char(27)//"[1;34m"//"===== FIXED SOURCE SOLVE ====="//char(27)//"[0m"
      src_ext = 0.0d0
      !Source RHS term
      do gg = 1,G
        do kk = 1, mesh%nz
          do jj = 1, mesh%ny
            do ii = 1, mesh%nx
              if (ii == ceiling(mesh%nx/2.0d0) .and. jj == jj .and. kk == kk) src_ext(1,ii + (jj-1)*mesh%nx + (kk-1)*mesh%nx*mesh%ny) = 1.0d0/(mesh%dx) !ceiling(mesh%nx/2.0d0)
            end do
          end do 
        end do
      end do

      do ii = 1, max_iter

        phi_prime = phi

        fission_rate = 0.0d0
        do gg = 1, G
          fission_rate(:) = fission_rate(:) + nu_Sigma_f(gg) * phi(gg,:)
        end do

        do gg = 1, G

          Sigma_s_L_sum = 0.0d0
          Sigma_s_U_sum = 0.0d0

          do jj = 1, G
            Sigma_s_L_sum(jj,:) = Sigma_s_L(jj,gg) * phi_prime(jj,:)
            Sigma_s_U_sum(jj,:) = Sigma_s_U(jj,gg) * phi(jj,:)
          end do

          S_f = (src_ext(gg,:) + sum(Sigma_s_L_sum, dim=1) + sum(Sigma_s_U_sum, dim=1) + chi(gg) * fission_rate(:)) * mesh%dV

          phi_ptr = phi(gg,:)
          call PCG_algorithm(A_all(gg)%aa, A_all(gg)%ja, A_all(gg)%ia, phi_ptr, S_f, PCG_mode, max_iter)
          phi(gg,:) = phi_ptr

        end do

        if (mod(ii, 10) == 0) then
            print *, ( sum(abs(phi(gg,:) - phi_prime(gg,:))), gg = 1, G )
        end if

        converged_inner =.true.
        do gg = 1, G
          if (sum(abs(phi(gg,:) - phi_prime(gg,:))) > 1.0d-6) then
            converged_inner =.false.
          end if
        end do
        if (converged_inner) exit
        
      end do

      call outpVTK_xyz_multi('fixed_flux', phi, mesh, G, use_adjoint)
      call outpVTK_xyz_multi('src_external', src_ext, mesh, G, use_adjoint)
      write(*,'(A)') char(27)//"[1;34m"//"===== FIXED SOLVER CONVERGED ====="//char(27)//"[0m"

    end if

    print *, ""
    if (use_adjoint) write(*,'(A)') char(27)//"[1;34m"//"===== ADJOINT EIGENVALUE SOLVE ====="//char(27)//"[0m"
    if (.not. use_adjoint) write(*,'(A)') char(27)//"[1;34m"//"===== PRIMAL EIGENVALUE SOLVE ====="//char(27)//"[0m"

    !------------------------------
    ! Outer power iteration with inner solve for each group
    !------------------------------
    outer: do ii = 2, max_iter
      S_f_iter = S_f
      inner: do
        do gg = 1, G

          ! Build scattering contributions into RHS:
          ! S_f = chi(gg) / K_eff * previous_fission_source + sum_downscatter(L*phi_prime) + sum_upscatter(U*phi)
          ! compute upscatter (i>j) contributions using phi_prime (previous inner iterate)
          Sigma_s_L_sum(:,:) = 0.0d0
          Sigma_s_U_sum(:,:) = 0.0d0

          do jj = 1, G
            ! Contributions for group gg are Sigma_s_L(jj,gg)*phi_prime(jj,:)
            Sigma_s_L_sum(jj,:) = Sigma_s_L(jj,gg) * phi_prime(jj,:)
            Sigma_s_U_sum(jj,:) = Sigma_s_U(jj,gg) * phi(jj,:)
            if (use_adjoint) Sigma_s_L_sum(jj,:) = Sigma_s_L(jj,gg) * phi(jj,:)
            if (use_adjoint) Sigma_s_U_sum(jj,:) = Sigma_s_U(jj,gg) * phi_prime(jj,:)
          end do

          S_f = ((chi(gg) / K_eff(ii-1)) * S_f_iter + sum(Sigma_s_L_sum, dim=1) + sum(Sigma_s_U_sum, dim=1))*mesh%dV
          if (use_adjoint) S_f = ((nu_sigma_f(gg)/K_eff(ii-1)) * S_f_iter  + sum(Sigma_s_L_sum,dim=1) + sum(Sigma_s_U_sum,dim=1))*mesh%dV
          
          phi_ptr = phi(gg, :)
          if (PCG_mode/=PETsc) call PCG_algorithm(A_all(gg)%aa, A_all(gg)%ja, A_all(gg)%ia, phi_ptr, S_f, PCG_mode, max_iter)
          if (PCG_mode==PETsc) call run_petsc(mesh%N, mesh%N, A_all(gg)%aa, A_all(gg)%ja, A_all(gg)%ia, phi_ptr, S_f)
          phi(gg, :) = phi_ptr

        end do

        converged_inner =.true.
        do gg = 1, G
          if (sum(abs(phi(gg,:) - phi_prime(gg,:))) >= 1.0d-8) then
            converged_inner =.false.
          end if
        end do
        if (converged_inner) exit inner
        phi_prime(:,:) = phi(:,:)

      end do inner

      ! compute new fission source and keff
      S_f = matmul(transpose(phi), nu_Sigma_f(1:G))
      if (use_adjoint) S_f = matmul(transpose(phi),chi(1:G))

      K_eff(ii) = K_eff(ii-1) * sum(S_f * mesh%dV) / sum(S_f_iter * mesh%dV)

      norm = sum(matmul(transpose(phi),nu_Sigma_f(1:G)) * mesh%dV)
      phi(:,:) = phi(:,:) / norm
        
      print *, "Iteration =", ii-1, "  Keff =", K_eff(ii-1)

      ! check outer convergence (relative keff change and fission source)
      if ( abs( (K_eff(ii) - K_eff(ii-1)) / K_eff(ii-1) ) < 1.0d-6 .and. maxval( abs( S_f - S_f_iter ) ) < 1.0d-6 ) exit
      !print*,S_f-S_f_iter
    end do outer

    if (use_adjoint) write(*,'(A)') char(27)//"[1;34m"//"===== ADJOINT CONVERGED ====="//char(27)//"[0m"
    if (.not. use_adjoint) write(*,'(A)') char(27)//"[1;34m"//"===== PRIMAL (FORWARD) CONVERGED ====="//char(27)//"[0m"

    phi = phi/(sum(matmul(transpose(phi), nu_Sigma_f(1:G)) *mesh%dV))  ! normalize flux

    call compute_cell_current(phi, D, mesh%nx, mesh%ny, mesh%nz, mesh%dx, mesh%dy, mesh%dz, Jx, Jy, Jz, G, use_adjoint)

    call compute_wall_currents(phi, mesh, D, BCs, Jx_min, Jx_max, Jy_min, Jy_max, Jz_min, Jz_max)

    do gg = 1, G
      call outpvtk_current('J_group'//trim(adjustl(itoa(gg)))//'.vtk', mesh, Jx(gg,:), Jy(gg,:), Jz(gg,:))
    end do

    print *, ""
    print '(A,3(F4.1,A))', "Cuboid Geometry: ", mesh%x_domain," x ",mesh%y_domain," x ",mesh%z_domain," cm"
    print '(A,3(I4,A),I3)', "     Mesh Cells: ", mesh%nx," x ",mesh%ny," x ",mesh%nz

    print *, ""
    print *, "=============================="
    print *, "Net Wall Currents: Groups 1-",trim(itoa(G))
    print *, "------------------------------"
    print "(A,*(E10.2))", "Jx_min:", ( sum(Jx_min(gg,:)), gg = 1, G )
    print "(A,*(E10.2))", "Jx_max:", ( sum(Jx_max(gg,:)), gg = 1, G )
    print "(A,*(E10.2))", "Jy_min:", ( sum(Jy_min(gg,:)), gg = 1, G )
    print "(A,*(E10.2))", "Jy_max:", ( sum(Jy_max(gg,:)), gg = 1, G )
    print "(A,*(E10.2))", "Jz_min:", ( sum(Jz_min(gg,:)), gg = 1, G )
    print "(A,*(E10.2))", "Jz_max:", ( sum(Jz_max(gg,:)), gg = 1, G )
    print *, "=============================="
    print *, ""

    ! Print outgoing group currents
    print *, "=============================="
    print *, "Outgoing Group Currents:"
    print *, "------------------------------"
    do gg = 1, G
      print "(A,I0,A,E12.4)", "Group ", gg, ": ", sum(Jx_min(gg,:)) + sum(Jx_max(gg,:)) + sum(Jy_min(gg,:)) + sum(Jy_max(gg,:)) + sum(Jz_min(gg,:)) + sum(Jz_max(gg,:))
    end do
    print *, "=============================="
    print *, ""

    ! Output VTK for visualization
    call outpVTK_xyz_multi('critical_flux',phi, mesh, G, use_adjoint)

    ! Compute production, absorption, and leakage
    P = sum(matmul(transpose(phi), nu_Sigma_f) * mesh%dV) / K_eff(ii)

    Abso = 0.0d0
    do gg = 1, G
        Abso = Abso + sum(Sigma_a(gg) * phi(gg,:) * mesh%dV)
    end do

    L = sum(Jx_min(:,:)) + sum(Jx_max(:,:)) + &
        sum(Jy_min(:,:)) + sum(Jy_max(:,:)) + &
        sum(Jz_min(:,:)) + sum(Jz_max(:,:))
    
    if (.not. use_adjoint) then
      print *,"=============================="
      print *, "Forward Balance Check"
      print *, "------------------------------"
      print *, "Production  =", P
      print *, "Absorption  =", Abso
      print *, "Leakage     =", L
      print *, "------------------------------"
      print *, "Balance     =", P - (Abso + L)
      print *, "=============================="
    end if

    print *, ""
    print *, "=============================="
    print *, "Keff =", K_eff(ii)
    print *, "=============================="

    print *, ""
    print *, "=============================="
    if (use_subcritical) call k_subcritical(mesh, MATs, K_eff(ii))
    print *, "=============================="

    print *, ""
    print *, "=============================="
    if (use_kinetics) call kinetic_solver(mesh, MATs, K_eff(ii))
    print *, "=============================="

    return
  end subroutine iter_xyz_multi

  pure function itoa(i) result(str)
    implicit none
    integer, intent(in) :: i
    character(len=16) :: str
    write(str,'(I0)') i
  end function itoa
  
end module

    ! print *, ""
    ! print *, "=============================="
    ! if (use_adjoint) print *, "Adjoint Material Properties"
    ! if (.not. use_adjoint) print *, "Primal Material Properties"
    ! print *, "------------------------------"
    ! print "(A)", "Sigma_s"
    ! do, ii=1,G
    !   write(*,"(100g15.5)") ( Sigma_s(ii,jj), jj=1,G)
    ! end do
    ! print "(A)", "Sigma a"
    ! print*,Sigma_a
    ! print "(A)", "nu Sigma f"
    ! print*, nu_Sigma_f
    ! print *, "=============================="