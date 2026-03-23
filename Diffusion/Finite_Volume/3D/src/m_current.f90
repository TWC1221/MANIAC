module m_current
  contains
  subroutine compute_cell_current(phi, Dg, nx, ny, nz, dx, dy, dz, Jx, Jy, Jz, G, use_adjoint)
    implicit none

    !-----------------------------
    ! Arguments
    !-----------------------------
    integer, intent(in) :: nx, ny, nz, G          ! mesh dimensions
    real(8), intent(in) :: dx, dy, dz          ! cell sizes

    logical, intent(in) :: use_adjoint

    ! Diffusion coefficients per group
    real(8), intent(in) :: Dg(:), phi(:,:)               ! size G

    ! Multigroup currents: J*(gg,row) => (G, N)
    real(8), intent(out) :: Jx(:,:), Jy(:,:), Jz(:,:)   ! size (G, N)

    !-----------------------------
    ! Locals
    !-----------------------------
    integer :: gg
    integer :: ii, jj, kk, row
    integer :: row_im1, row_ip1, row_jm1, row_jp1, row_km1, row_kp1
    real(8) :: dphidx, dphidy, dphidz
    real(8) :: Dloc

    !-----------------------------
    ! Main loops: over groups and cells
    !-----------------------------
    do gg = 1, G
      Dloc = Dg(gg)
      if (use_adjoint) Dloc = -Dloc

      do kk = 1, nz
        do jj = 1, ny
          do ii = 1, nx

            row = ii + (jj-1)*nx + (kk-1)*nx*ny

            !=========================
            ! dphi/dx for group gg
            !=========================
            if (ii == 1) then
                row_ip1 = (ii+1) + (jj-1)*nx + (kk-1)*nx*ny
                dphidx = (phi(gg, row_ip1) - phi(gg, row)) / dx
            else if (ii == nx) then
                row_im1 = (ii-1) + (jj-1)*nx + (kk-1)*nx*ny
                dphidx = (phi(gg, row) - phi(gg, row_im1)) / dx
            else
                row_im1 = (ii-1) + (jj-1)*nx + (kk-1)*nx*ny
                row_ip1 = (ii+1) + (jj-1)*nx + (kk-1)*nx*ny
                dphidx = (phi(gg, row_ip1) - phi(gg, row_im1)) / (2.0d0*dx)
            end if

            !=========================
            ! dphi/dy
            !=========================
            if (jj == 1) then
                row_jp1 = ii + (jj   )*nx + (kk-1)*nx*ny
                dphidy = (phi(gg, row_jp1) - phi(gg, row)) / dy
            else if (jj == ny) then
                row_jm1 = ii + (jj-2)*nx + (kk-1)*nx*ny
                dphidy = (phi(gg, row) - phi(gg, row_jm1)) / dy
            else
                row_jm1 = ii + (jj-2)*nx + (kk-1)*nx*ny
                row_jp1 = ii + (jj   )*nx + (kk-1)*nx*ny
                dphidy = (phi(gg, row_jp1) - phi(gg, row_jm1)) / (2.0d0*dy)
            end if

            !=========================
            ! dphi/dz
            !=========================
            if (kk == 1) then
                row_kp1 = ii + (jj-1)*nx + (kk   )*nx*ny
                dphidz = (phi(gg, row_kp1) - phi(gg, row)) / dz
            else if (kk == nz) then
                row_km1 = ii + (jj-1)*nx + (kk-2)*nx*ny
                dphidz = (phi(gg, row) - phi(gg, row_km1)) / dz
            else
                row_km1 = ii + (jj-1)*nx + (kk-2)*nx*ny
                row_kp1 = ii + (jj-1)*nx + (kk   )*nx*ny
                dphidz = (phi(gg, row_kp1) - phi(gg, row_km1)) / (2.0d0*dz)
            end if

            !print*,(phi(gg, row_kp1) - phi(gg, row_km1)) / (2.0d0*dz)

            !=========================
            ! J_g = -D_g * grad(phi_g)
            !=========================

            Jx(gg, row) = -Dloc * dphidx
            Jy(gg, row) = -Dloc * dphidy
            Jz(gg, row) = -Dloc * dphidz

          end do
        end do
      end do
    end do
  end subroutine compute_cell_current

end module